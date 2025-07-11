#include "../inc/util.hpp"
#include <fstream>

void write_csv(const std::string& filename, const std::vector<SimulationResult>& results) {
    std::ofstream file(filename);
    file << "NumberOfTubes,TubeHeight,RadiatorLength,TubeWidth,"
         << "FinDistance,FinHeight,FinWidth,RadiatorHeight,QValue\n";

    for (const auto& r : results) {
        file << r.NumberOfTubes << "," << r.TubeHeight << "," << r.RadiatorLength << ","
             << r.TubeWidth << "," << r.FinDistance << "," << r.FinHeight << ","
             << r.FinWidth << "," << r.RadiatorHeight << "," << r.QValue << "\n";
    }
}
