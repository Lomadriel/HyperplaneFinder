#ifndef HYPERPLANEFINDERD4_LINE_H
#define HYPERPLANEFINDERD4_LINE_H

#include <iostream>
#include <vector>

class Line {
public:
    Line(const std::vector<unsigned int>& points);

    Line(const std::initializer_list<unsigned int>& points);

    Line(std::vector<unsigned int>&& points);

    Line(const Line& line) = default;

    Line(Line&& line) = default;

    virtual ~Line() = default;

    std::vector<unsigned int>::size_type size() const;

    bool contains(const unsigned int point) const;

    bool isEmpty() const;

    bool isLine() const;

    bool isPoint() const;

    unsigned int getPoint(const unsigned int index) const;

    const std::vector<unsigned int>& getPoints() const;

    Line addLine(const Line& line) const;

    Line addScalar(const unsigned int scalar) const;

    Line intersects(const Line& line) const;

    Line unions(const Line& line) const;

    Line complementSymmetricDifference(const std::vector<unsigned int> set, const Line& line) const;

    bool isInclude(const Line& line) const;

    /* Operators */
    Line& operator=(const Line&) = default;

    bool operator==(const Line& rhs) const;

    bool operator!=(const Line& rhs) const;

    bool operator>(const Line& rhs) const;

    bool operator<(const Line& rhs) const;

    bool operator<=(const Line& rhs) const;

    bool operator>=(const Line& rhs) const;
private:
    std::vector<unsigned int> points;
};

std::ostream &operator<<(std::ostream &os, const Line& line);

#endif //HYPERPLANEFINDERD4_LINE_H
