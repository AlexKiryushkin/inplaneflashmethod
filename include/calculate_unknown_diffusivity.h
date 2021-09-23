
#include <vector>

#include "simulation_parameters.h"

namespace ipfm
{

IPFM_API double calculateUnknownDiffusivity(const SimulationParameters& params,
                                            std::vector<double> times,
                                            std::vector<double> temperatures);

} // namespace ipfm
