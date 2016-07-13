#include <iostream>
#include <algorithm>
#include <chrono>

#include <Line.h>
#include <PointsGeometry.h>

using namespace std;

int main() {
    std::vector<Line> lines;
    lines.push_back(Line({0, 1, 2}));
    lines.push_back(Line({3, 4, 5}));
    lines.push_back(Line({6, 7, 8}));
    lines.push_back(Line({0, 3, 6}));
    lines.push_back(Line({1, 4, 7}));
    lines.push_back(Line({2, 5, 8}));

    PointsGeometry geometry(9, std::move(lines));

    auto begin = chrono::high_resolution_clock::now();

    std::vector<Line> hyperplanes = geometry.findHyperplanes();
    cout << hyperplanes.size() << endl;
    std::sort(hyperplanes.begin(), hyperplanes.end());

    std::vector<Line> veldkampLines = geometry.findVeldkampLines(hyperplanes);
    cout << veldkampLines.size() << endl;

    std::vector<Line> hyperplanesX3 = geometry.computeHyperplanes();
    cout << hyperplanesX3.size() << endl;
    std::sort(hyperplanesX3.begin(), hyperplanesX3.end());

    auto end = chrono::high_resolution_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    cout << "Time: " << ms << endl;

    return 0;
}
