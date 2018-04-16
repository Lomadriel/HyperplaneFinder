#ifndef HYPERPLANEFINDER_EQUATIONSOLVER_HPP
#define HYPERPLANEFINDER_EQUATIONSOLVER_HPP

#include <bitset>
#include <tuple>
#include <vector>
#include "math_utility.hpp"

namespace segre {
	/*
	 * This type represents a solution of the equation.
	 * Each couple should be different from (0, 0).
	 * Moreover each coefficient should be in [0, 3[.
	 */
	template <std::size_t Dimension>
	using EquationSolution = std::array<std::array<unsigned char, 2>, Dimension>;

	template <std::size_t Dimension>
	std::vector<EquationSolution<Dimension>> resolveEquation(const std::array<unsigned char, math::pow(2, Dimension)>& coefficient);

	template <std::size_t Dimension>
	inline bool isSolution(const std::array<unsigned char, math::pow(2, Dimension)>& coefficient,
	                       const std::array<std::array<unsigned char, 2>, Dimension>& variables);

	template <std::size_t Dimension, std::size_t NumberOfPoints>
	std::bitset<NumberOfPoints> solutionsToHyperplane(const std::vector<EquationSolution<Dimension>>& solutions);

	template <std::size_t Dimension>
	std::size_t computeSolutionValue(const EquationSolution<Dimension>& solution);

	template <size_t Dimension>
	std::vector<EquationSolution<Dimension>> resolveEquation(const std::array<unsigned char, math::pow(2, Dimension)>& coefficient) {
		std::vector<EquationSolution<Dimension>> solutions{};

		std::array<std::array<unsigned char, 2>, Dimension> variables{};

		for (std::size_t i = 0; i < 3; ++i) {
			for (std::size_t j = 0; j < 3; ++j) {
			    for (std::size_t k = 0; k < 3; ++k) {
			        for (std::size_t l = 0; l < 3; ++l) {
			            for (std::size_t m = 0; m < 3; ++m) {
			                for (std::size_t n = 0; n < 3; ++n) {
			                    for (std::size_t o = 0; o < 3; ++o) {
			                        for (std::size_t p = 0; p < 3; ++p) {
			                        	variables[0][0] = i;
				                        variables[0][1] = j;
				                        variables[1][0] = k;
				                        variables[1][1] = l;
				                        variables[2][0] = m;
				                        variables[2][1] = n;
				                        variables[3][0] = o;
				                        variables[3][1] = p;

				                        if (variables[0][0] == 0 && variables[0][1] == 0) {
					                        continue;
				                        }

				                        if (variables[1][0] == 0 && variables[1][1] == 0) {
					                        continue;
				                        }

				                        if (variables[2][0] == 0 && variables[2][1] == 0) {
					                        continue;
				                        }

				                        if (variables[3][0] == 0 && variables[3][1] == 0) {
					                        continue;
				                        }
				                        if (isSolution(coefficient, variables)) {
				                        	solutions.push_back(variables);
				                        }
			                        }
			                    }
			                }
			            }
			        }
			    }
			}
		}

		return solutions;
	}

	template <size_t Dimension>
	inline bool isSolution(const std::array<unsigned char, math::pow(2, Dimension)>& coefficient,
	                const std::array<std::array<unsigned char, 2>, Dimension>& variables) {
		int value = 0;
		for (std::size_t i = 0; i < 2; ++i) {
			for (std::size_t j = 0; j < 2; ++j) {
				for (std::size_t k = 0; k < 2; ++k) {
					for (std::size_t l = 0; l < 2; ++l) {
						value += coefficient[l + k * 2 + j * 2 * 2 + i * 2 * 2 * 2] * variables[0][i] * variables[1][j] * variables[2][k] * variables[3][l];
					}
				}
			}
		}

		return value % 3 == 0;
	}

	template <std::size_t Dimension, std::size_t NumberOfPoints>
	std::bitset<NumberOfPoints> solutionsToHyperplane(const std::vector<EquationSolution<Dimension>>& solutions) {
		std::bitset<NumberOfPoints> hyperplane{};

		for (const EquationSolution<Dimension>& solution : solutions) {
			hyperplane[computeSolutionValue(solution)] = 1;
		}

		return hyperplane;
	}

	template <std::size_t Dimension>
	std::size_t computeSolutionValue(const EquationSolution<Dimension>& solution) {
		std::size_t value = 0;

		std::size_t i = 0;
		for (const std::array<unsigned char, 2>& couple : solution) {
			/*if (couple == {{1, 0}} || couple == {{2, 0}}) {
			}*/

			if (couple == std::array<unsigned char, 2>{{0, 1}} || couple == std::array<unsigned char, 2>{{0, 2}}) {
				value += 1 * math::pow(4, i);
			} else if (couple == std::array<unsigned char, 2>{{1, 1}} || couple == std::array<unsigned char, 2>{{2, 2}}) {
				value += 2 * math::pow(4, i);
			} else if (couple == std::array<unsigned char, 2>{{1, 2}} || couple == std::array<unsigned char, 2>{{2, 1}}) {
				value += 3 * math::pow(4, i);
			}

			++i;
		}

		return value;
	}
}

#endif // HYPERPLANEFINDER_EQUATIONSOLVER_HPP