
#pragma once

#include <cmath>
#include <ostream>

#include "import.h"

namespace ipfm
{

struct IPFM_API SimulationParameters
{
    double R;
    double rv;
    double ri;
    double ro;
    double alpha;
    double cp;
    double rho;
};

inline IPFM_API double defaultAlpha() { return 5.6; }

} // namespace ipfm

