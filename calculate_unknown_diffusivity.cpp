
#include "calculate_unknown_diffusivity.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

bool validateParams(const ipfm::SimulationParameters& params, std::string& errorString)
{
    if (params.cp <= 0.0) { errorString = "specific heat is incorrect"; return false; }
    if (params.rho <= 0.0) { errorString = "density is incorrect"; return false; }
    if (params.tmax <= 0.0) { errorString = "maximum time is incorrect"; return false; }
    if (params.alpha <= 0.0) { errorString = "heat exchange parameter is incorrect"; return false; }
    if (params.R <= 0.0) { errorString = "sample radius is incorrect"; return false; }
    if (params.rv <= 0.0) { errorString = "sensor radius is incorrect"; return false; }
    if (params.ri <= 0.0) { errorString = "inner heat radius is incorrect"; return false; }
    if (params.ro <= 0.0) { errorString = "outer heat radius is incorrect"; return false; }

    if (params.ro <= params.ri) { errorString = "outer heat radius should be greater than inner heat radius"; return false; }
    if (params.R < params.ro) { errorString = " sample radius should be greater than or equal to outer heat radius"; return false; }

    return true;
}

double J0(double x)
{
    return std::cyl_bessel_j(0, x);
}

double J1(double x)
{
    return std::cyl_bessel_j(1, x);
}

/**
 * @brief Finds root between leftX and rightX for bessel J1 function. Uses simple binary division for finding the root.
 * @return double 
 */
double findRoot(double leftX, double rightX, double leftFunctionValue, double rightFuntionValue)
{
    constexpr auto eps = 1e-7;

    const auto middleX = (leftX + rightX) / 2;
    const auto middleFunctionValue = J1(middleX);
    if (std::fabs(leftFunctionValue) < eps)
    {
        return leftX;
    }
    else if (std::fabs(rightFuntionValue) < eps)
    {
        return rightX;
    }
    else if (std::fabs(middleFunctionValue) < eps)
    {
        return middleX;
    }
    else if ( leftFunctionValue * middleFunctionValue < 0 )
    {
        return findRoot(leftX, middleX, leftFunctionValue, middleFunctionValue);
    }
    else if ( middleFunctionValue * rightFuntionValue < 0 )
    {
        return findRoot(middleX, rightX, middleFunctionValue, rightFuntionValue);
    }
    else
    {
        std::cout << "error\n";
        return 0.0;
    }
}

/**
 * @brief finds first @nRoots roots of J1 Bessel function
 */
std::vector<double> findFirstKRoots(std::size_t nRoots)
{
    std::vector<double> roots;

    const auto step = 0.01;
    double currentX = step;
    double currentFunctionValue = J1(currentX);
    for (std::size_t i = 1U; i <= nRoots; ++i)
    {
        auto nextX = currentX + step;
        double nextFunctionValue = J1(nextX);

        // go while function does not change its sign on the examined segment [currentX, nextX]
        while (currentFunctionValue * nextFunctionValue > 0)
        {
            currentX = nextX;
            currentFunctionValue = nextFunctionValue;

            nextX += step;
            nextFunctionValue = J1(nextX);
        }

        // add new root
        roots.push_back(findRoot(currentX, nextX, currentFunctionValue, nextFunctionValue));

        // go to next segment
        currentX = nextX;
        currentFunctionValue = nextFunctionValue;
    }

    return roots;
}

double phi(const ipfm::SimulationParameters& params, double betta)
{
    const auto rv = params.rv;
    const auto ro = params.ro;
    const auto ri = params.ri;
    const auto R = params.R;

    return J1(betta * rv / R) * ( ro * J1(betta * ro / R) - ri * J1(betta * ri / R) ) / 
        ( betta * betta * rv * J0(betta) * J0(betta) );
}

double V(const ipfm::SimulationParameters& params, const std::vector<double>& roots, double t, double lambda)
{
    constexpr double pi = 3.14159265358;
    const auto R = params.R;
    const auto ro = params.ro;
    const auto ri = params.ri;
    const auto rv = params.rv;
    const auto coeff = 4 * R * R / (ro * ro - ri * ri);
    const auto m = params.alpha / 2 / pi / R / R / params.cp / params.rho;

    double value{ 1.0 };
    for (const auto root : roots)
    {
        value += coeff * phi(params, root) * std::exp(-root * root * lambda * t / R / R);
    }
    value *= std::exp(-m * lambda * t / R / R);
    return value;
}

double calculateMaximumTemperatureTime(const ipfm::SimulationParameters& params, double lambda, double deltaT, double startT, double maxTime)
{
    const std::size_t nRoots = 50U;
    static auto roots = findFirstKRoots(nRoots);

    // find period where temperature starts to decrease
    double t = startT;
    double prevTemp = V(params, roots, t - deltaT, lambda);
    double currentTemp = V(params, roots, t, lambda);
    double nextTemp = V(params, roots, t + deltaT, lambda);
    while (t < maxTime)
    {
        if (currentTemp >= nextTemp && currentTemp >= prevTemp)
        {
            break;
        }

        t += deltaT;
        prevTemp = currentTemp; 
        currentTemp = nextTemp;
        nextTemp = V(params, roots, t + deltaT, lambda);
    }
    if (t >= maxTime)
    {
        return t;
    }

    // y = ax^2 + bx + c
    double a = (nextTemp - 2 * currentTemp + prevTemp) / 2 / deltaT / deltaT;
    double b = (nextTemp - prevTemp) / 2 / deltaT;
    double c = currentTemp;
    if (std::fabs(a) < std::numeric_limits<double>::epsilon())
    {
        return t;
    }

    double root = - b / 2 / a;
    return t  + root;
}

} // namespace


namespace ipfm
{

double calculateUnknownDiffusivity(const SimulationParameters& params)
{
    std::string error;
    if (!validateParams(params, error))
    {
        throw std::runtime_error{ error };
    }

    const auto timeSlices = 10000U;
    const auto deltaT = params.tmax / timeSlices;

    double minDiffusivity = 0.01;
    double maxDiffusivity = 3000;

    double minTime = calculateMaximumTemperatureTime(params, minDiffusivity, deltaT, deltaT, 1.1 * params.tmax);
    double maxTime = calculateMaximumTemperatureTime(params, maxDiffusivity, deltaT, deltaT, 1.1 * params.tmax);
    if (maxTime > params.tmax)
    {
        throw std::runtime_error{ "Provided time of reaching maximum temperature is too small. Probably there is a mistake in provided data." };
    }
    else if (minTime < params.tmax)
    {
        throw std::runtime_error{ "Provided time of reaching maximum temperature is too large. Probably there is a mistake in provided data." };
    }
    
    constexpr auto eps = 1e-7;
    double meanDiffusivity = 0.5 * (minDiffusivity + maxDiffusivity);
    double meanTime = calculateMaximumTemperatureTime(params, meanDiffusivity, deltaT, deltaT, 1.1 * params.tmax);
    while (std::fabs(meanTime - params.tmax) > eps)
    {
        if (meanTime > params.tmax)
        {
            minTime = meanTime;
            minDiffusivity = meanDiffusivity;
        }
        else
        {
            maxTime = meanTime;
            maxDiffusivity = meanDiffusivity;
        }
        meanDiffusivity = 0.5 * (minDiffusivity + maxDiffusivity);
        meanTime = calculateMaximumTemperatureTime(params, meanDiffusivity, deltaT, 5 * deltaT, 1.1 * params.tmax);
    }

    std::ofstream outFile{ "temp.csv" };
    outFile.imbue( std::locale( outFile.getloc(), new std::numpunct_byname<char>("de_DE.utf8") ) );
    const auto numIntervals = 1000U;
    const std::size_t nRoots = 30U;
    static auto roots = findFirstKRoots(nRoots);
    for (std::size_t idx{}; idx < numIntervals; ++idx)
    {
        const auto t = params.tmax * 4 * idx / numIntervals;
        outFile << t << ";" << V(params, roots, t, meanDiffusivity) << "\n";
    }

    return meanDiffusivity;
}

} // namespace ipfm

