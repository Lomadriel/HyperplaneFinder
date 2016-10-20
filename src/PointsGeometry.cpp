#include <vector>
#include <algorithm>
#include <numeric>
#include <thread>

#include <CombinationGenerator.h>
#include <PointsGeometry.h>
#include <future>

PointsGeometry::PointsGeometry(const unsigned int numberOfPoints,
                               const std::vector<Line>& lines)
    : elements(numberOfPoints)
    , lines(lines) {

    std::iota(elements.begin(), elements.end(), 0);
}

PointsGeometry::PointsGeometry(const unsigned int numberOfPoints,
                               std::vector<Line>&& lines)
        : elements(numberOfPoints)
        , lines(lines) {

    std::iota(elements.begin(), elements.end(), 0);
}

std::vector<Line> PointsGeometry::findHyperplanes() const {
    std::vector<Line> hyperplanes;

    CombinationGenerator gen;
    std::vector<unsigned int> currentCombination;

    for (unsigned int k = 0; k < elements.size(); ++k) {
        gen.initialize(static_cast<unsigned int>(elements.size()), k);

        while (!gen.isFinished()) {
            currentCombination = gen.nextCombination();
            Line line(currentCombination);

            if (isHyperplane(line)) {
                hyperplanes.push_back(std::move(line));
            }
        }
    }
    return hyperplanes;
}

bool PointsGeometry::isHyperplane(const Line& potentialHyperplane) const {
    for (const Line& line : lines) {
        Line intersect = potentialHyperplane.intersects(line);

        if (intersect.isEmpty()) {
            return false;
        }

        if (intersect.isLine()) {
            if (!potentialHyperplane.isInclude(line)) {
                return false;
            }
        }
    }

    return true;
}

std::vector<Line> PointsGeometry::computeCartesianProduct(const std::vector<Line>& geometry,
                                                          const unsigned int numberOfPoints,
                                                          const Line& line) {
    std::vector<Line> result(geometry);

    for (unsigned int i = 1; i < line.size(); ++i) {
        for (const Line& line1 : geometry) {
            result.push_back(line1.addScalar(i * numberOfPoints));
        }
    }

    std::vector<unsigned int> pts;
    for (unsigned int i = 0; i < numberOfPoints; ++i) {
        pts.reserve(line.size());
        for (unsigned int j = 0; j < line.size(); ++j) {
            pts[j] = i + j * numberOfPoints;
        }

        result.push_back(std::move(Line(std::move(pts))));
    }

    return result;
}

std::vector<Line> PointsGeometry::computeHyperplanes() const {
    std::vector<Line> hyperplanes = std::move(findHyperplanes());
    std::sort(hyperplanes.begin(), hyperplanes.end());

    return computeHyperplanes(std::move(hyperplanes));
}

std::vector<Line> PointsGeometry::computeHyperplanes(const std::vector<Line>& veldkampPoints) const {
    return computeHyperplanes(veldkampPoints, findVeldkampLines(veldkampPoints));
}

std::vector<Line> PointsGeometry::computeHyperplanes(const std::vector<Line>& veldkampPoints,
                                                     const std::vector<Line>& veldkampLines) const {
    std::vector<Line> hyperplanes;

    for (int i = static_cast<int>(veldkampLines.size()); i--;) {
        std::vector<Line> hypers = getHyperplanes(veldkampPoints, veldkampLines[i]);
        std::vector<std::vector<Line>> permutations = computePermutations(hypers);

        for (const std::vector<Line>& permutation : permutations) {
            Line hyperplane = permutation[0].addScalar(0 * static_cast<unsigned int>(elements.size()));
            for (int j = 1; j < permutation.size(); ++j) {
                Line shifted = permutation[j].addScalar(j * static_cast<unsigned int>(elements.size()));
                hyperplane = hyperplane.addLine(shifted);
            }
            hyperplanes.push_back(hyperplane);
        }
    }

    Line fullGeometry(elements);

    for (auto i = veldkampPoints.size(); i--;) {
        Line hyperplane = veldkampPoints[i];
        hyperplane.addLine(veldkampPoints[i].addScalar(static_cast<unsigned int>(elements.size())));
        hyperplane.addLine(fullGeometry.addScalar(2 * static_cast<unsigned int>(elements.size())));

        hyperplanes.push_back(hyperplane);

        hyperplane = veldkampPoints[i];
        hyperplane.addLine(fullGeometry.addScalar(static_cast<unsigned int>(elements.size())));
        hyperplane.addLine(veldkampPoints[i].addScalar(2 * static_cast<unsigned int>(elements.size())));

        hyperplanes.push_back(hyperplane);

        hyperplane = fullGeometry;
        hyperplane.addLine(veldkampPoints[i].addScalar(static_cast<unsigned int>(elements.size())));
        hyperplane.addLine(veldkampPoints[i].addScalar(2 * static_cast<unsigned int>(elements.size())));

        hyperplanes.push_back(hyperplane);
    }

    return hyperplanes;
}

std::vector<Line> PointsGeometry::computeHyperplanesD4(const std::vector<Line>& veldkampPoints,
                                                     const std::vector<Line>& veldkampLines) const {
    std::vector<Line> hyperplanes;

    for (int i = static_cast<int>(veldkampLines.size()); i--;) {
        std::vector<Line> hypers = getHyperplanes(veldkampPoints, veldkampLines[i]);
        std::vector<std::vector<Line>> permutations = computePermutations(hypers);

        for (const std::vector<Line>& permutation : permutations) {
            Line hyperplane = permutation[0].addScalar(0 * static_cast<unsigned int>(elements.size()));
            for (int j = 1; j < permutation.size(); ++j) {
                Line shifted = permutation[j].addScalar(j * static_cast<unsigned int>(elements.size()));
                hyperplane = hyperplane.addLine(shifted);
            }
            hyperplanes.push_back(hyperplane);
        }
    }

    Line fullGeometry(elements);


    for (auto i = veldkampPoints.size(); i--;) {
        Line hyperplane = veldkampPoints[i];
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(2 * static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(fullGeometry.addScalar(3 * static_cast<unsigned int>(elements.size())));

        hyperplanes.push_back(hyperplane);

        hyperplane = veldkampPoints[i];
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(fullGeometry.addScalar(2 * static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(3 * static_cast<unsigned int>(elements.size())));

        hyperplanes.push_back(hyperplane);

        hyperplane = veldkampPoints[i];
        hyperplane = hyperplane.addLine(fullGeometry.addScalar(static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(2 * static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(3 * static_cast<unsigned int>(elements.size())));

        hyperplanes.push_back(hyperplane);

        hyperplane = fullGeometry;
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(2 * static_cast<unsigned int>(elements.size())));
        hyperplane = hyperplane.addLine(veldkampPoints[i].addScalar(3 * static_cast<unsigned int>(elements.size())));

        hyperplanes.push_back(hyperplane);
    }

    return hyperplanes;
}

std::vector<Line> PointsGeometry::findVeldkampLines(const std::vector<Line>& veldkampPoints) const {
    std::vector<Line> veldkampLines;

    CombinationGenerator gen;
    gen.initialize(static_cast<unsigned int>(veldkampPoints.size()), 2);

    std::vector<unsigned int> currentCombination;

    while (!gen.isFinished()) {
        currentCombination = gen.nextCombination();
        const Line& h1 = veldkampPoints[currentCombination[0]];
        const Line& h2 = veldkampPoints[currentCombination[1]];
        const Line& hyperplane = h1.complementSymmetricDifference(elements, h2);

        unsigned int index3 =
                std::lower_bound(veldkampPoints.begin(), veldkampPoints.end(), hyperplane) - veldkampPoints.begin();

        if (currentCombination[1] < index3) {
            Line veldkampLine({currentCombination[0], currentCombination[1], index3});
            veldkampLines.push_back(veldkampLine);
        }
    }

    return veldkampLines;
}

std::pair<std::vector<Line>, std::vector<Line>> PointsGeometry::findVeldkampLinesP4(
        const std::vector<Line>& veldkampPoints) const {
    std::vector<Line> veldkampLineExc; // veldkamp lines supposed exceptional
    std::vector<Line> veldkampLineProj; // veldkamp lines projected

    CombinationGenerator gen;
    gen.initialize(static_cast<unsigned int>(veldkampPoints.size()), 2);

    std::vector<unsigned int> currentCombination;

    while (!gen.isFinished()) {
        currentCombination = gen.nextCombination();

        std::vector<std::pair<Line, unsigned int>> sameCore;
        const Line& h1 = veldkampPoints[currentCombination[0]];
        const Line& h2 = veldkampPoints[currentCombination[1]];
        const Line& intersection12 = h1.intersects(h2);

        for (unsigned int i = 0; i < veldkampPoints.size(); ++i) {
            const Line& intersection1i = h1.intersects(veldkampPoints[i]);
            const Line& intersection2i = h2.intersects(veldkampPoints[i]);

            if (intersection12 == intersection1i
                && intersection12 == intersection2i
                && intersection1i == intersection2i) {
                sameCore.push_back(std::make_pair(veldkampPoints[i], i));
            }
        }

        CombinationGenerator gen2;
        gen2.initialize(static_cast<unsigned int>(sameCore.size()), 2);
        std::vector<unsigned int> currentCombination2;

        while (!gen2.isFinished()) {
            currentCombination2 = gen2.nextCombination();

            const Line& ha = sameCore[currentCombination2[0]].first;
            const Line& hb = sameCore[currentCombination2[1]].first;
            const Line& intersectionAB = ha.intersects(hb);

            if (intersection12 == intersectionAB) {
                Line line = Line({currentCombination[0],
                                 currentCombination[1],
                                 sameCore[currentCombination2[0]].second,
                                 sameCore[currentCombination2[1]].second});

                if (std::find(veldkampLineExc.begin(), veldkampLineExc.end(), line) == veldkampLineExc.end()
                        && std::find(veldkampLineProj.begin(), veldkampLineProj.end(), line) == veldkampLineProj.end()) {
                    if (sameCore.size() == 2) {
                        veldkampLineProj.push_back(std::move(line));
                    } else {
                        veldkampLineExc.push_back(std::move(line));
                    }
                }
            }
        }
    }

    return std::make_pair(std::move(veldkampLineExc), std::move(veldkampLineProj));
}

void PointsGeometry::distinguishVeldkampLines(std::pair<std::vector<Line>, std::vector<Line>>& lines) const {
    // lines.first = exc lines
    // lines.second = proj lines
    size_t firstsz = lines.first.size();
    for (size_t i = 0; i < firstsz; ++i) {
        bool trulyExceptional = false;
        for (size_t j = 0, secondsz = lines.second.size(); j < secondsz; ++j) {
            if (lines.first[i].intersects(lines.second[j]).size() > 1) {
                trulyExceptional = true; // this exceptional line is truly exceptional.
                break;
            }
        }

        // still not truly exceptional.
        if (!trulyExceptional) { // this line is simply a projective one finally...
            // swap the exc line at the end to pop and put in proj.
            std::swap(lines.first[i--], lines.first[--firstsz]); // notice the post and pre decrementation.
            lines.second.push_back(lines.first.back());          
            lines.first.erase(lines.first.end() - 1, lines.first.end());
        }
    }
}

Line PointsGeometry::computeComplementHyperplane(const Line& h1,
                                                 const Line& h2,
                                                 const std::vector<unsigned int>& elements1) const {
    std::vector<unsigned int> complement;

    unsigned int i = 0, j = 0;
    for (auto& element : elements1) {
        if (i < h1.size() && h1.getPoint(i) == element) {
            ++i;
            if (j < h2.size() && h2.getPoint(j) == element) {
                ++j;
                complement.push_back(element);
            }
        } else if (j < h2.size() && h2.getPoint(j) == element) {
            ++j;
        } else {
            complement.push_back(element);
        }
    }

    return Line(std::move(complement));
}

std::vector<Line> PointsGeometry::getHyperplanes(const std::vector<Line>& veldkampPoints,
                                                 const Line& veldkampLine) {
    std::vector<Line> hyperplanes;

    for (unsigned int i = 0; i < veldkampLine.size(); ++i) {
        hyperplanes.push_back(veldkampPoints[veldkampLine.getPoint(i)]);
    }

    return hyperplanes;
}

PointsGeometry::Entry PointsGeometry::getEntry(const std::vector<Line> &lines,
                                               const Line &hyperplane) const {
    Entry entry;
    std::vector<Line> includedLines;
    for (const Line& line : lines) {
        if (hyperplane.isInclude(line)) {
            includedLines.push_back(line);
        }
    }

    entry.nbrPoints = static_cast<unsigned int>(hyperplane.size());
    entry.nbrLines = static_cast<unsigned int>(includedLines.size());

    if (entry.nbrLines == 0) {
        entry.pointsOfOrder[0] = entry.nbrPoints;
    } else {
        std::vector<unsigned int> integers;
        for (const Line& includedLine : includedLines) {
            integers.reserve(integers.size() + includedLine.size());
            std::vector<unsigned int> pts = includedLine.getPoints();
            integers.insert(integers.end(), pts.begin(),
                            pts.end());
        }

        std::sort(integers.begin(), integers.end());

        int count = 0;
        for (unsigned int i = 0, n = static_cast<unsigned int>(elements.size() * 3); i < n; ++i) {
            /*unsigned int countOfI =
                    static_cast<unsigned int>(std::count(integers.begin(), integers.end(), i));*/
            unsigned int countOfI = std::upper_bound(integers.begin(), integers.end(), i)
                                    - std::lower_bound(integers.begin(), integers.end(), i);

            if (countOfI > 0) {
                ++count;
                entry.pointsOfOrder[countOfI]++;
            }
        }

        int zero = entry.nbrPoints - count;
        if (zero != 0) {
            entry.pointsOfOrder[0] = zero;
        }
    }
    return entry;
}

std::vector<PointsGeometry::Entry> PointsGeometry::makeTable(const std::vector<Line>& veldkampPoints) const {
    std::vector<Entry> result;

    std::function<std::vector<Entry>(int)> compute = [&veldkampPoints, this](int a) -> std::vector<Entry> {
        std::vector<Entry> entries;

        for (int i = a; i < veldkampPoints.size(); i += 8) {
            PointsGeometry::Entry entry = this->getEntry(veldkampPoints, veldkampPoints[i]);

            std::vector<PointsGeometry::Entry>::iterator it = std::find(entries.begin(), entries.end(), entry);
            if (it == entries.end()) {
                entry.count = 1;
                entries.push_back(entry);
            } else {
                int index = it - entries.begin();
                entries[index].count += 1;
            }

            if (i > 0 && i % 1000 == 0) {
                std::cout << "Thread: " << a << ", " << i << " lines maked" << std::endl;
                break;
            }
        }
        return entries;
        };

    std::array<std::future<std::vector<Entry>>, 8> threadResult;
    for (int i = 0; i < threadResult.size(); ++i) {
        threadResult[i] = std::async(compute, i);
    }

    std::array<std::vector<Entry>, 8> finalResult;
    for (int i = 0; i < threadResult.size(); ++i) {
        finalResult[i] = threadResult[i].get();
    }

    return result;
}


std::vector<std::vector<Line>> PointsGeometry::computePermutations(const std::vector<Line>& veldkampLines) {
    std::vector<std::vector<Line>> permutations{veldkampLines};
    std::vector<Line> pickupList(veldkampLines);

    int* idx = new int[pickupList.size() + 1];
    for (int i = 0; i < pickupList.size() + 1; ++i) {
        idx[i] = i;
    }

    register int i = 1, j;
    while (i < pickupList.size()) {
        idx[i]--;
        j = i % 2 * idx[i];

        std::swap(pickupList[j], pickupList[i]);
        permutations.push_back(pickupList);

        i = 1;
        while (!idx[i]) {
            idx[i] = i;
            i++;
        }

    }

    delete[] idx;

    return permutations;
}

std::ostream& operator<<(std::ostream &os, const PointsGeometry::Entry& entry) {
    os << "Entry{" <<
    "Ps: " << entry.nbrPoints <<
    ", Ls: " << entry.nbrLines <<
    ", Order={";

    for(auto iterator = entry.pointsOfOrder.begin(); iterator != entry.pointsOfOrder.end();) {
        os << iterator->first << "=" << iterator->second;
        if (++iterator != entry.pointsOfOrder.end()) {
            os << ", ";
        }
    }

    os << "}";

    os << ", Crd: " << entry.count << '}';

    return os;
}
