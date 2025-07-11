#pragma once

#include <string>
#include <vector>

struct SimulationResult {
    double NumberOfTubes, TubeHeight, RadiatorLength, TubeWidth;
    double FinDistance, FinHeight, FinWidth;
    double RadiatorHeight, QValue;
};

void write_csv(const std::string& filename, const std::vector<SimulationResult>& results);
