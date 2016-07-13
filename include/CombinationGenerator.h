#ifndef HYPERPLANEFINDERD4_COMBINATIONGENERATOR_H
#define HYPERPLANEFINDERD4_COMBINATIONGENERATOR_H

#include <vector>

class CombinationGenerator {
public:
    CombinationGenerator() = default;
    virtual ~CombinationGenerator() = default;

    void initialize(const int n, const int k);

    void initialize(const int n, const int k, const unsigned long long numCombinations);

    std::vector<unsigned int> nextCombination();

    static unsigned long long computeBinomialCoefficient(int n, int k);

    unsigned long long getNumLeft() const;

    bool isFinished() const;

private:
    int n;
    int k;
    unsigned long long m_iNumCombinations;
    unsigned long long m_iNumLeft;
    std::vector<unsigned int> m_vecCurrCombination;

};

#endif //HYPERPLANEFINDERD4_COMBINATIONGENERATOR_H
