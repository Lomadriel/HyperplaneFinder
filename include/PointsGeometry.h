#ifndef HYPERPLANEFINDERD4_POINTSGEOMETRY_H
#define HYPERPLANEFINDERD4_POINTSGEOMETRY_H

#include <vector>
#include <map>
#include <Line.h>

class PointsGeometry {
public:
    struct Entry {
        template <typename Map>
        bool map_compare (Map const &lhs, Map const &rhs) const {
            return lhs.size() == rhs.size()
                   && std::equal(lhs.begin(), lhs.end(),
                                 rhs.begin());
        }

        bool operator==(const Entry& entry) const {
            return nbrPoints == entry.nbrPoints &&
                   nbrLines == entry.nbrLines &&
                   map_compare(pointsOfOrder, entry.pointsOfOrder);
        }

        friend std::ostream& operator<<(std::ostream &os, const PointsGeometry::Entry& entry);
    private:
        int nbrPoints;
        int nbrLines;
        std::map<int, int> pointsOfOrder;
        int count;

        friend class PointsGeometry;
    };

    typename PointsGeometry::Entry;

    PointsGeometry(const unsigned int numberOfPoints, const std::vector<Line>& lines);

    PointsGeometry(const unsigned int numberOfPoints, std::vector<Line>&& lines);

    std::vector<Line> findHyperplanes() const;

    static std::vector<Line> computeCartesianProduct(const std::vector<Line>& geometry,
                                                     const unsigned int numberOfPoints,
                                                     const Line& line);

    std::vector<Line> computeHyperplanes() const;

    std::vector<Line> computeHyperplanes(const std::vector<Line>& veldkampPoints) const;

    std::vector<Line> computeHyperplanes(const std::vector<Line>& veldkampPoints,
                                         const std::vector<Line>& veldkampLines) const;

    std::vector<Line> findVeldkampLines(const std::vector<Line>& veldkampPoints) const;

    std::pair<std::vector<Line>, std::vector<Line>> findVeldkampLinesD4(const std::vector<Line>& veldkampPoints) const;

    void distinguishVeldkampLines(std::pair<std::vector<Line>, std::vector<Line>>& lines) const;

private:
    bool isHyperplane(const Line& potentialHyperplane) const;

    Line computeComplementHyperplane(const Line& h1,
                                     const Line& h2,
                                     const std::vector<unsigned int>& element) const;

    static std::vector<Line> getHyperplanes(const std::vector<Line>& veldkampPoints,
                                            const Line& veldkampLine);

    Entry getEntry(const std::vector<Line>& lines, const Line& hyperplane) const;

    static std::vector<std::vector<Line>> computePermutations(const std::vector<Line>& veldkampLines);

private:
    std::vector<unsigned int> elements;
    std::vector<Line> lines;
};

std::ostream& operator<<(std::ostream &os, const PointsGeometry::Entry& entry);

#endif //HYPERPLANEFINDERD4_POINTSGEOMETRY_H
