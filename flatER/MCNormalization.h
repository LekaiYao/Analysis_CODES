#ifndef FLATER_MCNORMALIZATION_H
#define FLATER_MCNORMALIZATION_H

#include <algorithm>
#include <cctype>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace MCNormalization {

struct Context {
    std::string system;
    std::string tree;
    std::string particle;
    std::string promptness;
};

struct Entry {
    std::string system;
    std::string tree;
    std::string particle;
    std::string promptness;
    int pthat = -1;
    std::string pathPattern;
    double xsecPb = 0.;
    double filterEfficiency = 0.;
    long long nGenerated = 0;

    double Weight() const
    {
        return xsecPb * filterEfficiency / static_cast<double>(nGenerated);
    }
};

inline std::string Trim(const std::string &value)
{
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

inline std::string Lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

inline std::vector<std::string> SplitCsvLine(const std::string &line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(Trim(field));
    return fields;
}

inline Context BuildContext(const std::string &system,
                            const std::string &tree,
                            const std::string &particle,
                            const std::string &promptness)
{
    Context context;
    context.system = Lower(Trim(system));
    context.tree = Lower(Trim(tree));

    if (context.tree == "ntmix") {
        const std::string particleLower = Lower(particle);
        if (particleLower.find("psi2s") != std::string::npos) context.particle = "psi2s";
        else if (particleLower.find("x3872") != std::string::npos) context.particle = "x3872";
        else throw std::runtime_error("Unknown ntmix MC particle tag: '" + particle + "'");

        context.promptness = Lower(promptness).find("nonprompt") != std::string::npos
                           ? "nonprompt" : "prompt";
    } else if (context.tree == "ntkp") {
        context.particle = "bplus";
        context.promptness = "inclusive";
    } else if (context.tree == "ntphi") {
        context.particle = "bs";
        context.promptness = "inclusive";
    } else if (context.tree == "ntkstar") {
        context.particle = "bzero";
        context.promptness = "inclusive";
    } else {
        throw std::runtime_error("No MC normalization particle mapping for tree '" + tree + "'");
    }
    return context;
}

class Registry {
public:
    Registry() = default;

    explicit Registry(const std::string &fileName)
    {
        Load(fileName);
    }

    void Load(const std::string &fileName)
    {
        entries_.clear();
        std::ifstream input(fileName);
        if (!input.is_open()) {
            throw std::runtime_error("Cannot open MC normalization table: " + fileName);
        }

        std::string line;
        int lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;

            const auto fields = SplitCsvLine(line);
            if (Lower(fields.empty() ? "" : fields[0]) == "system") continue;
            if (fields.size() != 9) {
                throw std::runtime_error("Expected 9 columns in " + fileName + ":" +
                                         std::to_string(lineNumber));
            }

            Entry entry;
            try {
                entry.system = Lower(fields[0]);
                entry.tree = Lower(fields[1]);
                entry.particle = Lower(fields[2]);
                entry.promptness = Lower(fields[3]);
                entry.pthat = std::stoi(fields[4]);
                entry.pathPattern = Lower(fields[5]);
                entry.xsecPb = std::stod(fields[6]);
                entry.filterEfficiency = std::stod(fields[7]);
                entry.nGenerated = std::stoll(fields[8]);
            } catch (const std::exception &error) {
                throw std::runtime_error("Invalid value in " + fileName + ":" +
                                         std::to_string(lineNumber) + " (" + error.what() + ")");
            }

            if (entry.system.empty() || entry.tree.empty() || entry.particle.empty() ||
                entry.promptness.empty() || entry.pathPattern.empty() || entry.pthat < 0 ||
                entry.xsecPb <= 0. || entry.filterEfficiency <= 0. ||
                entry.filterEfficiency > 1. || entry.nGenerated <= 0) {
                throw std::runtime_error("Non-physical or empty value in " + fileName + ":" +
                                         std::to_string(lineNumber));
            }

            for (const auto &previous : entries_) {
                if (previous.system == entry.system && previous.tree == entry.tree &&
                    previous.particle == entry.particle &&
                    previous.promptness == entry.promptness &&
                    previous.pathPattern == entry.pathPattern) {
                    throw std::runtime_error("Duplicate MC normalization path pattern in " +
                                             fileName + ":" + std::to_string(lineNumber));
                }
            }
            entries_.push_back(entry);
        }

        if (entries_.empty()) {
            throw std::runtime_error("MC normalization table contains no entries: " + fileName);
        }
    }

    Entry Resolve(const std::string &fileName, const Context &context) const
    {
        const std::string pathLower = Lower(fileName);
        std::vector<const Entry *> matches;
        for (const auto &entry : entries_) {
            if (entry.system == context.system && entry.tree == context.tree &&
                entry.particle == context.particle &&
                entry.promptness == context.promptness &&
                pathLower.find(entry.pathPattern) != std::string::npos) {
                matches.push_back(&entry);
            }
        }

        if (matches.size() != 1) {
            std::ostringstream message;
            message << "Expected exactly one MC normalization row, found " << matches.size()
                    << " for file '" << fileName << "' (system=" << context.system
                    << ", tree=" << context.tree << ", particle=" << context.particle
                    << ", promptness=" << context.promptness << ")";
            throw std::runtime_error(message.str());
        }

        ValidatePthatInPath(pathLower, *matches.front());
        return *matches.front();
    }

private:
    static void ValidatePthatInPath(const std::string &path, const Entry &entry)
    {
        const std::regex pthatExpression("(pthat|phat)[-_]?([0-9]+)",
                                         std::regex_constants::icase);
        bool found = false;
        for (std::sregex_iterator match(path.begin(), path.end(), pthatExpression), end;
             match != end; ++match) {
            found = true;
            const int pthatInPath = std::stoi((*match)[2].str());
            if (pthatInPath != entry.pthat) {
                throw std::runtime_error("pThat mismatch: table has " +
                                         std::to_string(entry.pthat) + ", path contains " +
                                         std::to_string(pthatInPath) + " in '" + path + "'");
            }
        }
        if (!found) {
            throw std::runtime_error("No pThat/phat token found in MC path: " + path);
        }
    }

    std::vector<Entry> entries_;
};

} // namespace MCNormalization

#endif
