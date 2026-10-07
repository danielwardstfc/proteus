#pragma once

#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stdexcept>

using namespace mfem;

/// Simple text input parser for MFEM solvers.
/// Supports lines of the form: key  value  // optional comment
/// Ignores section headers and comment lines (starting with *, #, or //).
/// Throws an error if a requested parameter is missing.

class InputParser
{
private:
    std::unordered_map<std::string, std::string> params;

    static std::string Trim(const std::string &s)
    {
        size_t start = s.find_first_not_of(" \t\r\n");
        size_t end = s.find_last_not_of(" \t\r\n");
        return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
    }

public:
    InputParser(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file)
        {
            throw std::runtime_error("Error: Cannot open input file: " + filename);
        }

        std::string line;
        while (std::getline(file, line))
        {
            line = Trim(line);
            if (line.empty() || line[0] == '*' || line[0] == '#' || line.rfind("//", 0) == 0)
                continue;

            // Remove inline comments
            size_t commentPos = line.find("//");
            if (commentPos != std::string::npos)
                line = line.substr(0, commentPos);

            line = Trim(line);
            if (line.empty())
                continue;

            std::istringstream iss(line);
            std::string key, value;
            if (iss >> key >> value)
            {
                params[key] = value;
            }
        }
    }

    bool Has(const std::string &key) const
    {
        return params.find(key) != params.end();
    }

    int GetInt(const std::string &key) const
    {
        auto it = params.find(key);
        if (it == params.end())
        {
            throw std::runtime_error("Missing integer parameter: " + key);
        }
        return std::stoi(it->second);
    }

    mfem::real_t GetReal(const std::string &key) const
    {
        auto it = params.find(key);
        if (it == params.end())
        {
            throw std::runtime_error("Missing real_t parameter: " + key);
        }
        return static_cast<mfem::real_t>(std::stod(it->second));
    }

    std::string GetString(const std::string &key) const
    {
        auto it = params.find(key);
        if (it == params.end())
        {
            throw std::runtime_error("Missing string parameter: " + key);
        }
        return it->second;
    }
};