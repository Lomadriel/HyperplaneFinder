#include <algorithm>

#include <Line.h>

Line::Line(const std::vector<unsigned int>& points)
    : points(points) {
}

Line::Line(const std::initializer_list<unsigned int>& points)
    : points(points) {
}

Line::Line(std::vector<unsigned int>&& points)
    : points(std::move(points)) {
}


std::vector<unsigned int>::size_type Line::size() const {
    return points.size();
}

bool Line::contains(const unsigned int point) const {
    return std::find(points.begin(), points.end(), point) != points.end();
}

bool Line::isEmpty() const {
    return size() == 0;
}

bool Line::isLine() const {
    return size() > 1;
}

bool Line::isPoint() const {
    return !isLine();
}

unsigned int Line::getPoint(const unsigned int index) const {
    return points[index];
}

const std::vector<unsigned int>& Line::getPoints() const {
    return points;
}

Line Line::addLine(const Line& line) const {
    std::vector<unsigned int> pts;

    for (unsigned int i = 0; i < line.size(); ++i) {
        pts.push_back(line.getPoint(i));
    }

    return Line(std::move(pts));
}

Line Line::addScalar(const unsigned int scalar) const {
    std::vector<unsigned int> pts;
    pts.reserve(points.size());

    for (unsigned int i = 0; i < size(); ++i) {
        pts.push_back(getPoint(i) + scalar);
    }

    return Line(std::move(pts));
}

Line Line::intersects(const Line& line) const {
    std::vector<unsigned int> pts;

    unsigned int point;
    for (unsigned int i = 0; i < line.size(); ++i) {
        point = line.getPoint(i);
        if (contains(point)) {
            pts.push_back(point);
        }
    }

    return Line(std::move(pts));
}

Line Line::unions(const Line& line) const {
    std::vector<unsigned int> pts = points;

    unsigned int point;
    for (unsigned int i = 0; i < line.size(); ++i) {
        point = line.getPoint(i);
        if (std::find(pts.begin(), pts.end(), point) != pts.end()) {
            pts.push_back(point);
        }
    }

    return Line(std::move(pts));
}

bool Line::isInclude(const Line& line) const {
    for (int i = 0; i < line.points.size(); ++i) {
        if (!contains(line.points[i])) {
            return false;
        }
    }

    return true;
}

bool Line::operator==(const Line &rhs) const {
    return points == rhs.points;
}

bool Line::operator!=(const Line &rhs) const {
    return !operator==(rhs);
}

bool Line::operator>(const Line& rhs) const {
    auto len1 = size();
    auto len2 = rhs.size();
    auto lim = std::min(len1, len2);
    std::vector<unsigned int>::size_type k = 0;

    while (k < lim) {
        if (points[k] > rhs.points[k]) {
            return true;
        } else if (points[k] < rhs.points[k]) {
            return false;
        }

        ++k;
    }

    return (len1 - len2) > 0;
}

bool Line::operator<(const Line& rhs) const {
    return operator>(rhs);
};

bool Line::operator<=(const Line& rhs) const {
    return !operator>(rhs);
}

bool Line::operator>=(const Line& rhs) const {
    return !operator<(rhs);
}

std::ostream &operator<<(std::ostream &os, const Line& line) {
    os << "Line={";
    for (unsigned int i = 0; i < line.size(); ++i) {
        os << line.getPoint(i);
        if (i + 1 != line.size()) {
            os << ", ";
        }
    }

    os << "}";

    return os;
}
