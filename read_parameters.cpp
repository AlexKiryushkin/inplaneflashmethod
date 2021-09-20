
#include "read_parameters.h"

#include <fstream>
#include <iostream>
#include <string>

namespace
{
    
bool parseDouble(const std::string & str, double & value)
{
    try
    {
        value = std::stod(str);
        return true;
    }
    catch(const std::exception&)
    {
        return false;
    }
}

std::string getParametersFilePath()
{
    auto currentFile = std::string(__FILE__);

    auto pos = currentFile.rfind('\\');
    if (pos == std::string::npos)
    {
        pos = currentFile.rfind('/');
    }

    currentFile.erase(pos);
    return currentFile + "/parameters.txt";
}

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

} // namespace

namespace ipfm
{

bool readParameters(SimulationParameters& params)
{
    try
    {
        static const std::string pathToFile = getParametersFilePath();
        std::ifstream inputFile{ pathToFile };
        if (!inputFile)
        {
            std::cout << "Parameters file was not found. Fallback to interactive mode.\n";
            return false;
        }

        params.alpha = defaultAlpha();
        std::string parsedString;
        while(std::getline(inputFile, parsedString))
        {
            const auto spacePos = parsedString.find(" ");
            const std::string fieldStr{ parsedString.data(), spacePos };
            const auto valueStr = parsedString.substr(spacePos);
            const double value = std::stod( valueStr );
            if (fieldStr == "R")
            {
                params.R = value;
            }
            else if (fieldStr == "rv")
            {
                params.rv = value;
            }
            else if (fieldStr == "ri")
            {
                params.ri = value;
            }
            else if (fieldStr == "ro")
            {
                params.ro = value;
            }
            else if (fieldStr == "cp")
            {
                params.cp = value;
            }
            else if (fieldStr == "rho")
            {
                params.rho = value;
            }
            else if (fieldStr == "tmax")
            {
                params.tmax = value;
            }
            else if (fieldStr == "alpha")
            {
                params.alpha = value;
            }
            else
            {
                std::cout << "Parameters file has incorrect format. Fallback to interactive mode.\n";
                return false;
            }
        }

        std::string errorString;
        if (false == validateParams(params, errorString))
        {
            std::cout << errorString << "Fallback to interactive mode.\n";
            return false;
        }

        std::cout << "Parameters file at path: \"" + pathToFile + "\" was read correctly.\n";
        std::cout << "Sample radius: " << params.R << "\n";
        std::cout << "Sensor radius: " << params.rv << "\n";
        std::cout << "Inner heat radius: " << params.ri << "\n";
        std::cout << "Outer heat radius: " << params.ro << "\n";
        std::cout << "Sample specific heat: " << params.cp << "\n";
        std::cout << "Sample density: " << params.rho << "\n";
        std::cout << "Time of maximum temperature: " << params.tmax << "\n";
        return true;
    }
    catch(const std::exception&)
    {
        std::cout << "Unexpected error occurred. Fallback to interactive mode.\n";
        return false;
    }

}

void runInteractive(SimulationParameters& params)
{
    std::string parameterStr;

    std::cout << "Interactive mode. \n";
    std::cout << "Enter sample radius. R = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, params.R))
    {
        std::cout << "You did not enter a valid number. Enter sample radius. R = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter sensor radius. rv = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, params.rv))
    {
        std::cout << "You did not enter a valid number. Enter sensor radius. rv = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter inner heat radius. ri = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, params.ri))
    {
        std::cout << "You did not enter a valid number. Enter inner heat radius. ri = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter outer heat radius. ro = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, params.ro))
    {
        std::cout << "You did not enter a valid number. Enter outer heat radius. ro = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter sample specific heat. cp = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, params.cp))
    {
        std::cout << "You did not enter a valid number. Enter sample specific heat. cp = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter sample density. rho = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, params.rho))
    {
        std::cout << "You did not enter a valid number. Enter sample density. rho = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter time of maximum temperature. tmax = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, params.tmax))
    {
        std::cout << "You did not enter a valid number. Enter time of maximum temperature. tmax = ";
        std::cin >> parameterStr;
    }

    std::string errorString;
    if (false == validateParams(params, errorString))
    {
        throw std::runtime_error{ errorString };
    }
}

} // namespace ipfm
