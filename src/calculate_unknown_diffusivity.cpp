
#include "calculate_unknown_diffusivity.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

namespace
{

bool validateParams(const ipfm::SimulationParameters& params, std::string& errorString)
{
    if (params.cp <= 0.0) { errorString = "specific heat is incorrect"; return false; }
    if (params.rho <= 0.0) { errorString = "density is incorrect"; return false; }
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

    const auto divident = J1(betta * rv / R) * ( ro * J1(betta * ro / R) - ri * J1(betta * ri / R) );
    const auto divisor = betta * betta * rv * J0(betta) * J0(betta);
    const auto result = divident / divisor;
    return result;
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

void normalize(std::vector<double>& values)
{
    const auto [minIter, maxIter] = std::minmax_element(std::begin(values), std::end(values));
    const auto min = *minIter;
    const auto max = *maxIter;
    for (auto&& value : values)
    {
        value = (value - min) / (max - min);
    }
}

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

};

struct my_functor : Functor<double>
{
    my_functor(const ipfm::SimulationParameters& params,
               const std::vector<double> & times,
               const std::vector<double> & temperatures,
               const std::vector<double> & roots,
               std::vector<double> & calculatedTemperatures)
        : Functor<double>(1,1),
            params_{ params },
            times_(times),
            temperatures_(temperatures),
            roots_(roots),
            calculatedTemperatures_(calculatedTemperatures)
        {}
    int operator()(const Eigen::VectorXd &diffusivity, Eigen::VectorXd &residual) const
    {
        for (std::size_t idx{}; idx < times_.size(); ++idx)
        {
            calculatedTemperatures_.at(idx) = V(params_, roots_, times_.at(idx), diffusivity(0));
        }
        normalize(calculatedTemperatures_);
        residual(0U) = std::transform_reduce(std::begin(temperatures_),
                                             std::end(temperatures_),
                                             std::begin(calculatedTemperatures_), 0.0, std::plus{},
        [](auto lhs, auto rhs) { return (lhs - rhs) * (lhs - rhs); }) / double(temperatures_.size());
        std::cout << "residual: " << residual(0U) << "\n";
        return 0;
    }

    const ipfm::SimulationParameters& params_;
    const std::vector<double> & times_;
    const std::vector<double> & temperatures_;
    const std::vector<double> & roots_;
    std::vector<double> & calculatedTemperatures_;
};

} // namespace


namespace ipfm
{

double calculateUnknownDiffusivity(const SimulationParameters& params,
                                   std::vector<double> times,
                                   std::vector<double> temperatures)
{
    std::string error;
    if (!validateParams(params, error))
    {
        throw std::runtime_error{ error };
    }

    times.pop_back();
    temperatures.pop_back();
    normalize(temperatures);

    Eigen::VectorXd diffusivity(1);
    diffusivity(0) = 0.75;

    const std::size_t nRoots = 100U;
    auto roots = findFirstKRoots(nRoots);
    std::vector<double> calculatedTemperatures(temperatures.size());
    const my_functor functional{ params, times, temperatures, roots, calculatedTemperatures };
    Eigen::NumericalDiff<my_functor> numDiff(functional);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor>,double> lm(numDiff);
    lm.parameters.maxfev = 2000;
    lm.parameters.xtol = 1.0e-10;

    int ret = lm.minimize(diffusivity);
    std::cout << lm.iter << std::endl;
    std::cout << ret << std::endl;


    std::ofstream outFile{ "temp.csv" };
    outFile.imbue( std::locale( outFile.getloc(), new std::numpunct_byname<char>("de_DE.utf8") ) );
    for (std::size_t idx{}; idx < times.size(); ++idx)
    {
        outFile << times.at(idx) << ";" << calculatedTemperatures.at(idx) << ";" << temperatures.at(idx) << "\n";
    }

    return diffusivity(0);
}

} // namespace ipfm

