#include <iostream>
#include <algorithm>
#include <chrono>

#include <Line.h>
#include <PointsGeometry.h>

using namespace std;

int main() {
    std::vector<Line> lines;
    lines.push_back(Line({0, 1, 2, 3}));
    lines.push_back(Line({4, 5, 6, 7}));
    lines.push_back(Line({8, 9, 10, 11}));
    lines.push_back(Line({12, 13, 14, 15}));
    lines.push_back(Line({0, 4, 8, 12}));
    lines.push_back(Line({1, 5, 9, 13}));
    lines.push_back(Line({2, 6, 10, 14}));
    lines.push_back(Line({3, 7, 11, 15}));

    PointsGeometry geometry(16, std::move(lines));

    auto begin = chrono::high_resolution_clock::now();

    std::vector<Line> hyperplanes = geometry.findHyperplanes();
    cout << hyperplanes.size() << endl;
    std::sort(hyperplanes.begin(), hyperplanes.end());

    std::pair<std::vector<Line>, std::vector<Line>> veldkampLines = geometry.findVeldkampLinesP4(hyperplanes);
    sort(veldkampLines.second.begin(), veldkampLines.second.end());
    veldkampLines.second.erase(unique(veldkampLines.second.begin(), veldkampLines.second.end()), veldkampLines.second.end());

    sort(veldkampLines.first.begin(), veldkampLines.first.end() );
    veldkampLines.first.erase(unique(veldkampLines.first.begin(), veldkampLines.first.end()), veldkampLines.first.end());

    cout << veldkampLines.second.size() << " "
         << veldkampLines.first.size() << " "
         << veldkampLines.second.size() + veldkampLines.first.size() << endl;
    geometry.distinguishVeldkampLines(veldkampLines);
    cout << veldkampLines.second.size() + veldkampLines.first.size() << endl;
    cout << "Second: " << veldkampLines.second.size() << endl;

    /*std::vector<Line> hyperplanesX3 = geometry.computeHyperplanes();
    cout << hyperplanesX3.size() << endl;
    std::sort(hyperplanesX3.begin(), hyperplanesX3.end());*/

    auto end = chrono::high_resolution_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    cout << "Time: " << ms << endl;

    return 0;
}
