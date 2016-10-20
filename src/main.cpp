#include <iostream>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iterator>

#include <Line.h>
#include <PointsGeometry.h>

/*
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

    unsigned int numberOfPoints = 16;
    PointsGeometry geometry(numberOfPoints, std::move(lines));

    auto begin = std::chrono::high_resolution_clock::now();

    std::string fileName("hyperplanes-D2.txt");
    std::ofstream output(fileName);

    std::vector<Line> hyperplanes = geometry.findHyperplanes();
    output << "Hyperplanes: " << hyperplanes.size() << std::endl;
    std::copy(hyperplanes.begin(), hyperplanes.end(),
              std::ostream_iterator<Line>(output, "\n"));

    std::pair<std::vector<Line>, std::vector<Line>> veldkampLines;

    for (int i = 0, j = 2; i < 1; ++i, ++j) {
        veldkampLines = geometry.findVeldkampLinesP4(hyperplanes);
        geometry.distinguishVeldkampLines(veldkampLines);
        output << "Number of exc: " << veldkampLines.first.size()
               << std::endl;
        std::copy(veldkampLines.first.begin(), veldkampLines.first.end(),
                  std::ostream_iterator<Line>(output, "\n"));

        output << "Number of proj: " << veldkampLines.second.size()
               << std::endl;
        std::copy(veldkampLines.second.begin(), veldkampLines.second.end(),
                  std::ostream_iterator<Line>(output, "\n"));

        output.close();
        fileName = "hyperplanes-D";
        fileName += std::to_string(j+1) + ".txt";
        output.open(fileName, std::ios_base::out);

        hyperplanes = geometry.computeHyperplanesD4(hyperplanes, veldkampLines.second);
        output << "\n" << "Hyperplanes: " << hyperplanes.size() << std::endl;
        std::copy(hyperplanes.begin(), hyperplanes.end(),
                  std::ostream_iterator<Line>(output, "\n"));

        lines = PointsGeometry::computeCartesianProduct(lines, numberOfPoints, Line({0, 1, 2, 3}));
        numberOfPoints *= 4;
        geometry = PointsGeometry(numberOfPoints, lines);
    }

    begin = std::chrono::high_resolution_clock::now();

    output.close();
    auto entries = geometry.makeTable(hyperplanes);
    output.open("hyperplane-table-D4.txt");
    std::copy(entries.begin(), entries.end(), std::ostream_iterator<PointsGeometry::Entry>(output, "\n"));

    auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << "Time: " << ms << std::endl;

    return 0;
}
*/

int main() {
    std::ifstream input("/home/lomadriel/Code/ClionProjects/HyperplaneFinderD4/build/hyperplanes-D4.txt");
    std::string s;
    std::getline(input, s, ':');
    std::getline(input, s, '\n');
    unsigned long max = std::stoul(s);
    std::vector<Line> lines;

    std::cout << "Loading " << max << " lines" << std::endl;
    Line line;
    int i = 0;
    while (i < max && input >> line) {
        lines.push_back(line);
        ++i;
        line = Line();

        if (i % 1000000 == 0) {
            std::cout << "Loaded " << i << "lines" << std::endl;
        }
    }

    input.close();

    std::vector<Line> geometryLines;
    geometryLines.push_back(Line({0, 1, 2, 3}));
    geometryLines.push_back(Line({4, 5, 6, 7}));
    geometryLines.push_back(Line({8, 9, 10, 11}));
    geometryLines.push_back(Line({12, 13, 14, 15}));
    geometryLines.push_back(Line({0, 4, 8, 12}));
    geometryLines.push_back(Line({1, 5, 9, 13}));
    geometryLines.push_back(Line({2, 6, 10, 14}));
    geometryLines.push_back(Line({3, 7, 11, 15}));

    geometryLines = PointsGeometry::computeCartesianProduct(geometryLines, 16*4, Line({0, 1, 2, 3}));
    PointsGeometry geometry(16*4, geometryLines);

    std::cout << "Making table" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<PointsGeometry::Entry> entries = geometry.makeTable(lines);
    auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << "Time: " << ms << std::endl;

    std::ofstream ofstream("output.txt");
    for (auto& entry : entries) {
        ofstream << entry << std::endl;
    }

    ofstream.close();
}
