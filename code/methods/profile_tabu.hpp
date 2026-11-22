#pragma once

#include <chrono>

#include "../representation/Solution.hpp"

/**
 * @brief Weight-profile TabuCol descent.
 * The colour classes of the (feasible) solution define a profile of class levels (the max
 * weight of each class). Lowering the level of one class gives a cheaper profile; whether a
 * colouring fits it is a list-colouring decision problem (vertex v may only go in classes of
 * level >= w(v)), solved by TabuCol conflict minimisation. Every success is a strictly
 * cheaper legal colouring; the phase ends when every one-class lowering failed or at the
 * deadline.
 */
void profile_tabu(Solution &solution,
                  const std::chrono::high_resolution_clock::time_point &deadline);
