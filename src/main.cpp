
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
    if (!fs::exists(csvPath))
    {
        throw std::runtime_error{ "File path at " + csvPath.generic_string() + " does not exists" };
    }

    std::ifstream csvFile{ csvPath.generic_string() };
    std::cout << "Start parsing csv file..." << std::endl;

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
    std::cout << "Parsing finished." << std::endl;
    std::cout << "Number of expeimental data points is" << timeMoments.size() << std::endl;
    std::cout << "=====================================" << std::endl;

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

fs::path parseInputFile(int argc, char** argv)
{
    if (argc != 2)
    {
        throw std::runtime_error{ "Incorrect number of arguments. Run the program like: ipfm --input-file=path_to_file" };
    }

    fs::path exePath{ argv[0] };
    std::string lineArg{ argv[1] };

    const auto foundIdx = lineArg.find('=');
    std::string argName = lineArg.substr(0U, foundIdx);
    if (argName != "--input-file")
    {
        throw std::runtime_error{ "Incorrect argument name. Run the program like: ipfm --input-file=path_to_file" };
    }

    std::string argValue = lineArg.substr(foundIdx + 1U);
    if (fs::path filePath{argValue}; filePath.is_absolute())
    {
        if (!fs::exists(filePath))
        {
            throw std::runtime_error{ "There is no file at path" + filePath.generic_string() };
        }

        return argValue;
    }

    fs::path filePath = exePath.parent_path() / argValue;
    if (!fs::exists(filePath))
    {
        throw std::runtime_error{ "There is no file at path" + filePath.generic_string() };
    }

    return filePath;
}

int main(int argc, char** argv)
{
    try
    {
        auto inputFile = parseInputFile(argc, argv);
        std::cout << "Reading experimental data from " << inputFile.generic_string() << " file" << std::endl;

        auto params = ipfm::SimulationParameters{ };
        if (!ipfm::readParameters(params))
        {
            ipfm::runInteractive(params);
        }

        auto [time, temperature] = parseCsv(inputFile.generic_string());
        sparse(time, 10000U);
        sparse(temperature, 10000U);
        std::cout << "For performance reasons leave only " << time.size() << " experimental points." << std::endl;

        const auto sz = time.size();
        std::vector<double> smoothedTemperature(sz), tmp1(sz), tmp2(sz);
        CppLowess::TemplatedLowess<std::vector<double>, double> dlowess;

        std::cout << "Start smoothing experimental data..." << std::endl;
        dlowess.lowess(time, temperature, 0.05, 0, 0.0, smoothedTemperature, tmp1, tmp2);
        std::cout << "Finished smoothing" << std::endl;

        sparse(smoothedTemperature, 1000U);
        sparse(time, 1000U);
        std::cout << "For performance reasons leave only " << time.size() << " smoothed points." << std::endl;

        const auto unknownDiffusivity = ipfm::calculateUnknownDiffusivity(params, time, smoothedTemperature);
        std::cout << "=================================================================\n" << std::endl;
        std::cout << "Unknown thermal diffusivity is: " << unknownDiffusivity << " [cm^2/s]" << std::endl;
        std::cout << "Unknown thermal conductivity is: " << unknownDiffusivity * params.cp * params.rho << " [W/(cm K)]" << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}