#ifndef HYPERPLANEFINDER_COMBINATIONGENERATOR_HPP
#define HYPERPLANEFINDER_COMBINATIONGENERATOR_HPP

#include <vector>
#include <cassert>
#include <algorithm>
#include <numeric>

class CombinationGenerator {
public:
	CombinationGenerator()
			: m_k(0)
			, m_n(0)
			, m_iNumCombinations(0)
			, m_iNumLeft(0)
			, m_vecCurrCombination{} {
	}

	void initialize(const unsigned int n, const unsigned int k) {
		initialize(n, k, computeBinomialCoefficient(n, k));
	}

	void initialize(const unsigned int n, const unsigned int k, const unsigned long long numCombinations) noexcept {
		assert(n >= k && n > 1);
		m_n = n;
		m_k = k;

		m_iNumLeft = m_iNumCombinations = numCombinations;

		m_vecCurrCombination.resize(static_cast<unsigned long>(k));
		std::iota(m_vecCurrCombination.begin(), m_vecCurrCombination.end(), 0);
	}

	void reset() noexcept {
		m_iNumLeft = m_iNumCombinations;
		std::sort(m_vecCurrCombination.begin(), m_vecCurrCombination.end());
	}

	const std::vector<unsigned int>& nextCombination() noexcept {
		if(m_iNumLeft == m_iNumCombinations) {
			--m_iNumLeft;
			return m_vecCurrCombination;
		}

		size_t i = m_k - 1;
		while(m_vecCurrCombination[i] == m_n - m_k + i) {
			--i;
		}

		++(m_vecCurrCombination[i]);

		for(size_t j = i; j < m_k; ++j) {
			m_vecCurrCombination[j] = static_cast<unsigned int>(m_vecCurrCombination[i] + j - i);
		}

		--m_iNumLeft;

		return m_vecCurrCombination;
	}

	static unsigned long long computeBinomialCoefficient(unsigned int n, unsigned int k) {
		assert(n >= k && n > 1);

		if(k > n - k) {
			k = n - k;
		}

		unsigned long long c = 1;

		for(size_t i = 1; i < k + 1; i++) {
			c *= n - (k - i);
			c /= i;
		}

		return c;
	}

	unsigned long long getNumLeft() const {
		return m_iNumLeft;
	}

	bool isFinished() const {
		return m_iNumLeft == 0;
	}

private:
	unsigned int m_k;
	unsigned int m_n;
	unsigned long long m_iNumCombinations;
	unsigned long long m_iNumLeft;
	std::vector<unsigned int> m_vecCurrCombination;

};

#endif //HYPERPLANEFINDER_COMBINATIONGENERATOR_HPP
