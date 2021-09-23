
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "calculate_unknown_diffusivity.h"
#include "read_parameters.h"

#include "lowess.h"

namespace fs = std::filesystem;

std::pair<std::vector<double>, std::vector<double>> parseCsv(const fs::path& csvPath)
{
    auto currentFilePath = fs::path(__FILE__);
    currentFilePath = currentFilePath.parent_path().parent_path();
    currentFilePath.append(csvPath.generic_string());

    if (!fs::exists(currentFilePath))
    {
        throw std::runtime_error{ "File path at " + currentFilePath.generic_string() + " does not exists" };
    }

    std::ifstream csvFile{ currentFilePath.generic_string() };

    // skip first line
    std::string dummy;
    std::getline(csvFile, dummy);

    std::vector<double> timeMoments;
    std::vector<double> temperatureValues;
    while(csvFile)
    {
        char dummyC{};
        double time{};
        double temperature{};
        csvFile >> time >> dummyC >> temperature;

        timeMoments.push_back(time);
        temperatureValues.push_back(temperature);
    }

    return { std::move(timeMoments), std::move(temperatureValues) };
}

void sparse(std::vector<double>& values, std::size_t maxSize)
{
    const auto step = values.size() / maxSize;
    if (step < 2)
    {
        return;
    }

    std::vector<double> newValues;
    newValues.reserve(maxSize);

    for (std::size_t idx{}; idx < values.size(); idx += step)
    {
        newValues.push_back(values.at(idx));
    }

    values = std::move(newValues);
}

int main(int argc, char** argv)
{
    try
    {
        auto params = ipfm::SimulationParameters{ };
        if (!ipfm::readParameters(params))
        {
            ipfm::runInteractive(params);
        }

        auto [time, temperature] = parseCsv("data.csv");
        sparse(time, 20000U);
        sparse(temperature, 20000U);

        const auto sz = time.size();
        std::vector<double> smoothedTemperature(sz), tmp1(sz), tmp2(sz);
        CppLowess::TemplatedLowess<std::vector<double>, double> dlowess;
        dlowess.lowess(time, temperature, 0.05, 0, 0.0, smoothedTemperature, tmp1, tmp2);
        sparse(smoothedTemperature, 2000U);
        sparse(time, 2000U);

        const auto unknownDiffusivity = ipfm::calculateUnknownDiffusivity(params, time, smoothedTemperature);
        std::cout << "Unknown thermal diffusivity is: " << unknownDiffusivity << std::endl;
        std::cout << "Unknown thermal conductivity is: " << unknownDiffusivity * params.cp * params.rho << std::endl;
        
        std::ofstream out{ "smoothed_data.csv" };
        out.imbue( std::locale( out.getloc(), new std::numpunct_byname<char>("de_DE.utf8") ) );
        for (std::size_t i{}; i < time.size() - 1U; ++i)
        {
            out << time.at(i) << ";" << smoothedTemperature.at(i) << "\n";
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}