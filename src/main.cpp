#include "../inc/radiator.hpp"
#include "../inc/util.hpp"
#include <vector>

int main() {
    std::vector<SimulationResult> results;
    const int iterations = 100;
    const double radLenMin = 0.0, radLenMax = 1.0;
    const double step = (radLenMax - radLenMin) / iterations;

    const double numTubes     = 60;
    const double tubeHeight   = 0.00156;
    const double tubeWidth    = 0.0246;
    const double finDistance  = 0.00158;
    const double finHeight    = 0.1188;
    const double finWidth     = 0.0246;

    for (int i = 0; i <= iterations; ++i) {
        double radLength = radLenMin + i * step;

        Radiator sim(numTubes, tubeHeight, radLength, tubeWidth,
                     finDistance, finHeight, finWidth);

        double q = sim.run();

        results.push_back({
            numTubes, tubeHeight, radLength, tubeWidth,
            finDistance, finHeight, finWidth,
            sim.getRadiatorHeight(), q
        });
    }

    write_csv("results.csv", results);
    return 0;
}
