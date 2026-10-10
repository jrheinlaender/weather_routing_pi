#include "supercpn/weather_routing/Providers.h"
#include "supercpn/weather_routing/Types.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace supercpn::weather_routing {
namespace {
double rowSpeed(const PolarRow& row, double angle) {
  angle = std::clamp(std::abs(angle), 0.0, 180.0);
  auto upper =
      std::lower_bound(row.GetPoints().begin(), row.GetPoints().end(), angle,
                       [](const PolarPoint& point, double value) {
                         return point.trueWindAngleDegrees < value;
                       });
  if (upper == row.GetPoints().begin()) return upper->boatSpeedKnots;
  if (upper == row.GetPoints().end())
    return row.GetPoints().back().boatSpeedKnots;
  const auto& high = *upper;
  const auto& low = *(upper - 1);
  const double span = high.trueWindAngleDegrees - low.trueWindAngleDegrees;
  const double factor =
      span <= 0.0 ? 0.0 : (angle - low.trueWindAngleDegrees) / span;
  return low.boatSpeedKnots +
         factor * (high.boatSpeedKnots - low.boatSpeedKnots);
}

double waveFactor(const std::vector<WavePerformancePoint>& response,
                  double height, double period) {
  if (response.empty()) return 1.0;
  std::vector<std::pair<double, double>> nearest;
  nearest.reserve(response.size());
  for (const auto& point : response) {
    const double distance =
        std::hypot((height - point.significantHeightMetres) / 2.0,
                   (period - point.periodSeconds) / 6.0);
    if (distance < 1e-9) return point.speedFactor;
    nearest.emplace_back(distance, point.speedFactor);
  }
  std::stable_sort(nearest.begin(), nearest.end());
  double weighted = 0.0;
  double weights = 0.0;
  for (std::size_t i = 0; i < std::min<std::size_t>(4, nearest.size()); ++i) {
    const double weight = 1.0 / (nearest[i].first * nearest[i].first);
    weighted += nearest[i].second * weight;
    weights += weight;
  }
  return std::clamp(weighted / weights, 0.0, 1.5);
}
}  // namespace

PolarPerformanceModel::PolarPerformanceModel(VesselConfiguration configuration)
    : configuration_(std::move(configuration)) {}

bool PolarPerformanceModel::validateProfile(const PerformanceProfile& profile,
                                            std::string* reason) {
  auto fail = [&](const std::string& message) {
    if (reason) *reason = message;
    return false;
  };
  if (profile.identity.empty()) return fail("profile identity is empty");
  if (profile.rows.empty()) return fail("profile has no wind-speed rows");
  double previousWind = -1.0;
  for (const auto& row : profile.rows) {
    if (!std::isfinite(row.GetTws()) || row.GetTws() <= previousWind)
      return fail("true wind speeds must be finite and strictly increasing");
    if (row.GetPoints().size() < 2)
      return fail("each polar row needs at least two angles");
    double previousAngle = -1.0;
    for (const auto& point : row.GetPoints()) {
      if (!std::isfinite(point.trueWindAngleDegrees) ||
          !std::isfinite(point.boatSpeedKnots) ||
          point.trueWindAngleDegrees < 0.0 ||
          point.trueWindAngleDegrees > 180.0 ||
          point.trueWindAngleDegrees <= previousAngle ||
          point.boatSpeedKnots < 0.0)
        return fail(
            "polar angles/speeds are invalid or not strictly increasing");
      previousAngle = point.trueWindAngleDegrees;
    }
    previousWind = row.GetTws();
  }
  if (!std::isfinite(profile.efficiency) || profile.efficiency <= 0.0)
    return fail("profile efficiency must be positive");
  for (const auto& point : profile.wavePerformance)
    if (!std::isfinite(point.significantHeightMetres) ||
        point.significantHeightMetres < 0.0 ||
        !std::isfinite(point.periodSeconds) || point.periodSeconds <= 0.0 ||
        !std::isfinite(point.speedFactor) || point.speedFactor < 0.0 ||
        point.speedFactor > 1.5)
      return fail("wave performance points are invalid");
  return true;
}

bool PolarPerformanceModel::valid(std::string* reason) const {
  const auto* sailing = profile(ProfileRole::SailOnly);
  if (!sailing) {
    if (reason) *reason = "a sail-only profile is required";
    return false;
  }
  for (const auto& item : configuration_.profiles)
    if (!validateProfile(item, reason)) return false;
  if (configuration_.propulsion.allowMotor &&
      configuration_.propulsion.configuredMotorSpeedKnots <= 0.0 &&
      !profile(ProfileRole::MotorOnly)) {
    if (reason)
      *reason = "motor mode needs a motor profile or configured motor speed";
    return false;
  }
  return true;
}

const PerformanceProfile* PolarPerformanceModel::profile(
    ProfileRole role, const std::string& identity) const {
  for (const auto& item : configuration_.profiles)
    if (item.role == role && (identity.empty() || item.identity == identity))
      return &item;
  return nullptr;
}

double PolarPerformanceModel::interpolate(const PerformanceProfile& profile,
                                          double tws, double twa) {
  if (profile.rows.empty()) return std::numeric_limits<double>::quiet_NaN();
  auto upper = std::lower_bound(
      profile.rows.begin(), profile.rows.end(), tws,
      [](const PolarRow& row, double value) { return row.GetTws() < value; });
  if (upper == profile.rows.begin()) return rowSpeed(*upper, twa);
  if (upper == profile.rows.end()) return rowSpeed(profile.rows.back(), twa);
  const auto& high = *upper;
  const auto& low = *(upper - 1);
  const double factor = (tws - low.GetTws()) / (high.GetTws() - low.GetTws());
  return rowSpeed(low, twa) +
         factor * (rowSpeed(high, twa) - rowSpeed(low, twa));
}

PerformanceCandidate PolarPerformanceModel::evaluateProfile(
    const PerformanceProfile& profileValue, double tws, double twa,
    const WaveSample& waves) const {
  double speed = interpolate(profileValue, tws, twa) * profileValue.efficiency;
  const double directionalEfficiency = twa <= 90.0
                                           ? configuration_.upwindEfficiency
                                           : configuration_.downwindEfficiency;
  speed *= directionalEfficiency;
  if (waves.available) {
    if (!profileValue.wavePerformance.empty())
      speed *= waveFactor(profileValue.wavePerformance,
                          waves.significantHeightMetres, waves.periodSeconds);
    else if (configuration_.wavePenaltyPerMetre > 0.0)
      speed *= std::max(0.0, 1.0 - waves.significantHeightMetres *
                                       configuration_.wavePenaltyPerMetre);
  }
  PerformanceCandidate result;
  result.valid = std::isfinite(speed) && speed > 0.0;
  result.role = profileValue.role;
  result.profileIdentity = profileValue.identity;
  result.sailPlan =
      static_cast<int>(&profileValue - configuration_.profiles.data());
  result.speedThroughWaterKnots = std::max(0.0, speed);
  result.mode = profileValue.role == ProfileRole::SailOnly
                    ? PropulsionMode::Sail
                : profileValue.role == ProfileRole::MotorSailing
                    ? PropulsionMode::MotorSail
                    : PropulsionMode::Motor;
  if (result.mode != PropulsionMode::Sail &&
      configuration_.propulsion.fuelConsumptionLitresPerHour)
    result.fuelLitresPerHour =
        *configuration_.propulsion.fuelConsumptionLitresPerHour;
  return result;
}

PerformanceCandidate PolarPerformanceModel::evaluate(
    PropulsionMode mode, ProfileRole role, const std::string& identity,
    double tws, double twa, const WaveSample& waves) const {
  if (const auto* selected = profile(role, identity)) {
    auto result = evaluateProfile(*selected, tws, twa, waves);
    result.mode = mode;
    return result;
  }
  PerformanceCandidate result;
  result.mode = mode;
  result.role = role;
  result.profileIdentity =
      identity.empty() ? "configured-constant-motor" : identity;
  if (mode == PropulsionMode::Motor && configuration_.propulsion.allowMotor &&
      configuration_.propulsion.configuredMotorSpeedKnots > 0.0) {
    result.valid = true;
    result.speedThroughWaterKnots =
        configuration_.propulsion.configuredMotorSpeedKnots;
  } else if (mode == PropulsionMode::MotorSail &&
             configuration_.propulsion.allowMotorSailing &&
             configuration_.propulsion.motorSailingBoostKnots > 0.0) {
    if (const auto* sail = profile(ProfileRole::SailOnly)) {
      result = evaluateProfile(*sail, tws, twa, waves);
      result.mode = PropulsionMode::MotorSail;
      result.role = ProfileRole::MotorSailing;
      result.profileIdentity = "sail-plus-configured-motor-assist";
      result.speedThroughWaterKnots +=
          configuration_.propulsion.motorSailingBoostKnots;
    }
  }
  if (result.valid && mode != PropulsionMode::Sail &&
      configuration_.propulsion.fuelConsumptionLitresPerHour)
    result.fuelLitresPerHour =
        *configuration_.propulsion.fuelConsumptionLitresPerHour;
  return result;
}

namespace {
double deg2rad(const double angle) {
  return angle / 180.0 * std::numbers::pi;
}

double rad2deg(const double angle) {
  return angle / std::numbers::pi * 180.0;
}
using OptimalAngles::Angle::STARBOARD_UPWIND;
using OptimalAngles::Angle::STARBOARD_DOWNWIND;
using OptimalAngles::Angle::PORT_DOWNWIND;
using OptimalAngles::Angle::PORT_UPWIND;

double iterSolveOptimalAngle(const PolarPoint& low, const PolarPoint& high,
                             const double twa_start) {
  double d_twa = deg2rad(high.trueWindAngleDegrees - low.trueWindAngleDegrees);
  double d_bs = high.boatSpeedKnots - low.boatSpeedKnots;
  double twa_low = deg2rad(low.trueWindAngleDegrees);
  double twa_next = twa_start;
  double twa;
  unsigned it = 0;
  do {
    twa = twa_next;
    /**
     * The range of atan() is -90 to 90 degrees, therefore true wind angles
     * must be corrected. Note that using atan2() does not help. Since d_bs can
     * be positive or negative, while the denominator is always positive, the
     * result of atan2() will always be in the first or second quadrant.
     */
    twa_next = atan(d_bs / ((twa - twa_low) * d_bs + d_twa * low.boatSpeedKnots));
    if (twa > std::numbers::pi / 2)
      twa_next += std::numbers::pi;
    if (twa > 3 * std::numbers::pi / 2)
      twa_next += std::numbers::pi;
  } while (std::fabs(twa_next - twa) > 1E-5 && ++it < 20);

  return twa_next;
}
}

const std::shared_ptr<OptimalAngles> PolarRow::GetOptimalAngles() const {
  if (m_optimal_angles->valid())
    return m_optimal_angles;
  if (m_points.empty())
    return nullptr;

  std::array<std::pair<std::vector<PolarPoint>::const_iterator, double>, 4> p_max =
    {std::pair{m_points.end(), 0.0}, std::pair{m_points.end(), 0.0},
     std::pair{m_points.end(), 0.0}, std::pair{m_points.end(), 0.0}};

  for (auto p = m_points.begin(); p != m_points.end(); ++p) {
    double twa = p->trueWindAngleDegrees;
    // Ensure that 180° is in section 1, to handle symmetric polars
    auto sector = OptimalAngles::Angle((twa - 1E-6) / 90);
    double upwind = (sector == STARBOARD_UPWIND ||
                     sector == PORT_UPWIND) ? 1 : -1;

    double vmg = upwind * cos(deg2rad(twa)) * p->boatSpeedKnots;
    if (vmg > p_max[size_t(sector)].second)
      p_max[size_t(sector)] = {p, vmg};
  }

  // interpolate the best vmg (as it often lies between entries in the polar)
  for (auto a : OptimalAngles::All) {
    if (m_points.back().trueWindAngleDegrees <= 180.0 && a > STARBOARD_DOWNWIND)
      break; // Symmetric polar

    const auto& p_mid = p_max[size_t(a)].first;
    if (p_mid == m_points.end())
      continue;

    double upwind = (a == STARBOARD_UPWIND || a == PORT_UPWIND) ? 1.0 : -1.0;
    double twa_l = deg2rad(p_mid->trueWindAngleDegrees);
    double twa_r = twa_l;
    auto p_l = (p_mid == m_points.begin()) ? p_mid : p_mid - 1;
    auto p_r = (p_mid == m_points.end() - 1) ? p_mid : p_mid + 1;
    double twa = twa_l;

    /**
     * Check gradient to decide whether extremum is left or right of max_twa
     * Note that it is possible to have an extremum on either side
     * The SailingVMG method does not cover this case
     */
    double d_bs_l  = p_mid->boatSpeedKnots - p_l->boatSpeedKnots;
    double d_bs_r  = p_r->boatSpeedKnots - p_mid->boatSpeedKnots;
    double d_twa_l = deg2rad(p_mid->trueWindAngleDegrees - p_l->trueWindAngleDegrees);
    double d_twa_r = deg2rad(p_r->trueWindAngleDegrees - p_mid->trueWindAngleDegrees);
    double vmg_l   = 0.0;
    double vmg_r   = 0.0;

    if (d_twa_l > 1E-6) {
      double grad_l = d_bs_l / d_twa_l * cos(twa) -
                      sin(twa) * p_mid->boatSpeedKnots;
      if (grad_l * upwind < 0.0)
        twa_l = iterSolveOptimalAngle(*p_l, *p_mid, twa - 0.25 * d_twa_l);
      vmg_l = upwind * cos(twa_l) * (p_mid->boatSpeedKnots - (twa - twa_l) / d_twa_l * d_bs_l);
    }
    if (d_twa_r > 1E-6) {
      double grad_r = d_bs_r / d_twa_r * cos(twa) -
                      sin(twa) * p_mid->boatSpeedKnots;
      if (grad_r * upwind > 0.0)
        twa_r = iterSolveOptimalAngle(*p_mid, *p_r, twa + 0.25 * d_twa_r);
      vmg_r = upwind * cos(twa_r) * (p_mid->boatSpeedKnots + (twa_r - twa) / d_twa_r * d_bs_r);
    }

    m_optimal_angles->at(a) = rad2deg(vmg_l > vmg_r ? twa_l : twa_r);
  }

  // For symmetric polars. And it assumes the points are sorted ascending
  if (m_points.back().trueWindAngleDegrees <= 180.0) {
    m_optimal_angles->at(PORT_UPWIND) =
      360.0 - m_optimal_angles->at(STARBOARD_UPWIND);
    m_optimal_angles->at(PORT_DOWNWIND) =
      360.0 - m_optimal_angles->at(STARBOARD_DOWNWIND);
  }

  return m_optimal_angles;
}

std::shared_ptr<OptimalAngles> PolarPerformanceModel::GetOptimalAngles(
    PropulsionMode mode, ProfileRole role, const std::string& identity,
    double tws, const WaveSample&) const {
  if (mode != PropulsionMode::Sail || role != ProfileRole::SailOnly)
      return {};

  const PerformanceProfile* prof = profile(role, identity);
  if (prof == nullptr)
      return {};
  if (prof->rows.empty())
    return {};

  // Search polar rows between which the wind speed lies
  const auto after = std::lower_bound(
      prof->rows.begin(), prof->rows.end(), tws,
      [](const PolarRow& row, double value) {
          return row.GetTws() < value;
      });

  if (after == prof->rows.begin())
    return after->GetOptimalAngles();
  if (after == prof->rows.end())
    return prof->rows.back().GetOptimalAngles();
  const auto before = after - 1;
  // No optimal angles can be found at zero wind speed
  if (before->GetTws() < 1E-9)
    return after->GetOptimalAngles();

  const auto opt_before = before->GetOptimalAngles();
  const auto opt_after = after->GetOptimalAngles();

  if (!opt_before->valid() || !opt_after->valid())
      return {};

  // Linear interpolation, could be improved by using 1st derivatives
  std::shared_ptr<OptimalAngles> result = std::make_shared<OptimalAngles>();
  double tws_ratio = (tws - before->GetTws()) /
                     (after->GetTws() - before->GetTws());

  for (const auto a : OptimalAngles::All)
      result->at(a) = std::lerp(opt_before->at(a), opt_after->at(a), tws_ratio);

  return result;
}

std::vector<PerformanceCandidate> PolarPerformanceModel::candidates(
    double tws, double twa, const WaveSample& waves,
    PropulsionMode previousMode, Duration previousModeDuration) const {
  std::vector<PerformanceCandidate> result;
  PerformanceCandidate sail;
  if (configuration_.propulsion.allowSailing) {
    for (const auto& sailing : configuration_.profiles)
      if (sailing.role == ProfileRole::SailOnly) {
        auto candidate = evaluateProfile(sailing, tws, twa, waves);
        if (!candidate.valid) continue;
        if (!sail.valid || candidate.speedThroughWaterKnots >
                               sail.speedThroughWaterKnots)
          sail = candidate;
        result.push_back(std::move(candidate));
      }
  }
  const bool minimumRunBlocksSail =
      previousMode != PropulsionMode::Sail &&
      previousModeDuration < configuration_.propulsion.minimumMotorRun;
  if (minimumRunBlocksSail)
    result.erase(std::remove_if(result.begin(), result.end(),
                                [](const auto& item) {
                                  return item.mode == PropulsionMode::Sail;
                                }),
                 result.end());

  const double threshold =
      configuration_.propulsion.motorBelowSailingSpeedKnots;
  const bool motorNeeded =
      minimumRunBlocksSail || !sail.valid || threshold <= 0.0 ||
      sail.speedThroughWaterKnots <
          threshold +
              (previousMode == PropulsionMode::Sail
                   ? 0.0
                   : configuration_.propulsion.crossoverHysteresisKnots);
  if (motorNeeded && configuration_.propulsion.allowMotorSailing) {
    auto motorSail = evaluate(PropulsionMode::MotorSail,
                              ProfileRole::MotorSailing, {}, tws, twa, waves);
    if (motorSail.valid) result.push_back(std::move(motorSail));
  }
  if (motorNeeded && configuration_.propulsion.allowMotor) {
    auto motor = evaluate(PropulsionMode::Motor, ProfileRole::MotorOnly, {},
                          tws, twa, waves);
    if (motor.valid) result.push_back(std::move(motor));
  }
  return result;
}

}  // namespace supercpn::weather_routing
