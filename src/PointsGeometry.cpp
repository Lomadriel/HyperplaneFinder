#include <vector>
#include <algorithm>
#include <numeric>

#include <CombinationGenerator.h>
#include <PointsGeometry.h>

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
    pts.reserve(line.size());
    for (unsigned int i = 0; i < numberOfPoints; ++i) {
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

std::pair<std::vector<Line>, std::vector<Line>> PointsGeometry::findVeldkampLinesD4(const std::vector<Line>& veldkampPoints) const {
    std::vector<Line> veldkampLineExc; // veldkamp lines supposed exceptional
    std::vector<Line> veldkampLineProj; // veldkamp lines projected

    CombinationGenerator gen;
    gen.initialize(static_cast<unsigned int>(veldkampPoints.size()), 2);

    std::vector<unsigned int> currentCombination;

    while (!gen.isFinished()) {
        currentCombination = gen.nextCombination();

        std::vector<Line> core;
        const Line& h1 = veldkampPoints[currentCombination[0]];
        const Line& h2 = veldkampPoints[currentCombination[1]];
        const Line& intersection12 = h1.intersects(h2);

        for (const Line& veldkampPoint : veldkampPoints) {
            const Line& intersection10 = h1.intersects(veldkampPoint);
            const Line& intersection20 = h2.intersects(veldkampPoint);

            if (intersection12 == intersection10 && intersection12 == intersection20 && intersection10 == intersection20) {
                core.push_back(veldkampPoint);
            }
        }

        CombinationGenerator gen2;
        gen2.initialize(static_cast<unsigned int>(veldkampPoints.size()), 2);
        std::vector<unsigned int> currentCombination2;

        while (!gen.isFinished()) {
            currentCombination2 = gen.nextCombination();

            const Line& ha = veldkampPoints[currentCombination2[0]];
            const Line& hb = veldkampPoints[currentCombination2[1]];
            const Line& intersectionAB = ha.intersects(hb);

            if (intersection12 == intersectionAB) {
                if (core.size() == 2) {
                    veldkampLineProj.push_back(std::move(Line({currentCombination[0],
                                                               currentCombination[1],
                                                               currentCombination2[0],
                                                               currentCombination2[1]})));
                } else {
                    veldkampLineExc.push_back(std::move(Line({currentCombination[0],
                                                               currentCombination[1],
                                                               currentCombination2[0],
                                                               currentCombination2[1]})));
                }
            }
        }
    }

    return std::make_pair(std::move(veldkampLineExc), std::move(veldkampLineProj));
};

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

std::vector<std::vector<Line>> PointsGeometry::computePermutations(const std::vector<Line>& veldkampLines) {
    std::vector<std::vector<Line>> permutations{veldkampLines};
    std::vector<Line> pickupList(veldkampLines);

    int idx[pickupList.size() + 1];
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
