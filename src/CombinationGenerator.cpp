#include <vector>
#include <cassert>
#include <CombinationGenerator.h>

void CombinationGenerator::initialize(const int n, const int k) {
    initialize(n, k, computeBinomialCoefficient(n, k));
}

void CombinationGenerator::initialize(const int n, const int k, const unsigned long long numCombinations) {
    assert(n >= k && n > 1);
    this->n = n;
    this->k = k;

    m_iNumLeft = m_iNumCombinations = numCombinations;

    m_vecCurrCombination.clear();

    for (unsigned int i = 0; i < k; i++) {
        m_vecCurrCombination.push_back(i);
    }
}

unsigned long long CombinationGenerator::computeBinomialCoefficient(int n, int k) {
    assert(n > k && n > 1);

    if (k > n - k) {
        k = n - k;
    }

    unsigned long long c = 1;

    for (int i = 1; i < k+1; i++) {
        c *= n - (k - i);
        c /= i;
    }

    return c;
}

unsigned long long CombinationGenerator::getNumLeft() const {
    return m_iNumLeft;
}

bool CombinationGenerator::isFinished() const {
    return m_iNumLeft == 0;
}

std::vector<unsigned int> CombinationGenerator::nextCombination() {
    if (m_iNumLeft == m_iNumCombinations) {
        m_iNumLeft -= 1;
        return m_vecCurrCombination;
    }

    int i = k - 1;
    while (m_vecCurrCombination[i] == n - k + i) {
        i--;
    }

    m_vecCurrCombination[i] += 1;

    for (int j = i; j < k; j++) {
        m_vecCurrCombination[j] = m_vecCurrCombination[i] + j - i;
    }

    m_iNumLeft -= 1;

    return m_vecCurrCombination;
}
