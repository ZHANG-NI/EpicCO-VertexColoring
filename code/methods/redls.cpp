#include "redls.hpp"

#include <algorithm>
#include <cassert>

#include "../utils/random_generator.hpp"
#include "../utils/utils.hpp"
#include "profile_tabu.hpp"

namespace {
// Anti-cycling attribute tabu (ls_atabu > 0): when a vertex leaves a colour in the
// conflict phase, it may not re-enter that colour until move counter g_moves passes
// g_tabu_until[vertex * (n + 1) + colour]. Tenure = ls_atabu x |conflicting vertices| +
// rand(0..9) moves (TabuCol rule). Aspiration: a tabu move that makes the colouring
// feasible is allowed.
std::vector<long> g_tabu_until;
long g_moves = 0;

inline bool atabu_on() { return Parameters::p->ls_atabu > 0; }

inline bool atabu(const int vertex, const int color) {
    return g_tabu_until[static_cast<size_t>(vertex) * (Graph::g->nb_vertices + 1) +
                        static_cast<size_t>(color)] > g_moves;
}

inline void atabu_set(const int vertex, const int old_color, const int nb_conflicting) {
    ++g_moves;
    if (old_color < 0) {
        return;
    }
    std::uniform_int_distribution<int> dis(0, 9);
    const long tenure =
        static_cast<long>(Parameters::p->ls_atabu * static_cast<double>(nb_conflicting)) +
        dis(rd::generator);
    g_tabu_until[static_cast<size_t>(vertex) * (Graph::g->nb_vertices + 1) +
                 static_cast<size_t>(old_color)] = g_moves + tenure;
}
} // namespace

void redls(Solution &best_solution, const bool verbose) {
    auto max_time = std::chrono::high_resolution_clock::now() +
                    std::chrono::seconds(Parameters::p->max_time_local_search);

    int64_t best_time = 0;

    if (atabu_on()) {
        const size_t sz = static_cast<size_t>(Graph::g->nb_vertices) *
                          static_cast<size_t>(Graph::g->nb_vertices + 1);
        if (g_tabu_until.size() != sz) {
            g_tabu_until.assign(sz, 0);
            g_moves = 0;
        }
        // entries of earlier calls expire
        g_moves += 1L << 30;
    }
    ProxiSolutionRedLS solution(best_solution);
    std::vector<bool> tabu_list(Graph::g->nb_vertices, false);
    long turn = 0;
    // slack injection (ls_inject > 0, offer S1_v7): after a window of turns without a new
    // local best (doubling after every injection without progress), apply the reduce step
    // to the current infeasible colouring so stuck conflicting vertices get slack
    const long inject_base =
        Parameters::p->ls_inject > 0
            ? std::max(10000L, static_cast<long>(Parameters::p->ls_inject *
                                                 static_cast<double>(Graph::g->nb_vertices)))
            : 0;
    long inject_window = inject_base;
    long last_progress = 0;
    while (not Parameters::p->time_limit_reached_sub_method(max_time) and
           turn < Parameters::p->nb_iter_local_search and
           best_solution.score_wvcp() != Parameters::p->target) {
        ++turn;
        if (solution.penalty() == 0) {
            while (improve_conflicts(solution, false, tabu_list)) {
                assert(solution.check_solution());
            }

            assert(solution.penalty() == 0);
            if (solution.score_wvcp() < best_solution.score_wvcp()) {
                last_progress = turn;
                inject_window = inject_base;
                best_solution = solution.solution();
                if (verbose) {
                    best_time = Parameters::p->elapsed_time(
                        std::chrono::high_resolution_clock::now());
                    print_result_ls(best_time, best_solution, turn);
                }
                if (solution.score_wvcp() == Parameters::p->target) {
                    if (verbose) {
                        print_result_ls(best_time, best_solution, turn);
                    }
                    return;
                }
            } else if (solution.score_wvcp() == best_solution.score_wvcp()) {
                best_solution = solution.solution();
            }
            std::fill(tabu_list.begin(), tabu_list.end(), false);
            move_heaviest_vertices(solution);
            assert(solution.check_solution());
            if (solution.penalty() == 0) {
                // for the case where moving the heaviest vertices to
                // a new color does not bring conflicts
                continue;
            }
        }

        if (inject_base > 0 and turn - last_progress > inject_window) {
            last_progress = turn;
            inject_window *= 2;
            move_heaviest_vertices(solution);
            if (solution.penalty() == 0) {
                continue;
            }
        }
        if (not improve_conflicts_and_score(
                solution, best_solution.score_wvcp(), tabu_list)) {

            while (improve_conflicts(solution, true, tabu_list)) {
                assert(solution.check_solution());
            }

            if (not solve_one_conflict_preserve_score(
                    solution, best_solution.score_wvcp(), tabu_list)) {
                if (not ejection_move(solution, best_solution.score_wvcp(), tabu_list)) {
                    // Increments edge weight
                    solution.increment_edge_weights();
                    if (not solution.conflict_edges().empty()) {
                        solve_one_conflict(solution, best_solution.score_wvcp(), tabu_list);
                    }
                }
            }
            assert(solution.check_solution());
        }
    }
    // gain-bounded gating (offer S2_v5): PT runs only when the largest PT gain seen so far
    // in the run could still bring this round below the incumbent (first calls always run)
    static int pt_calls = 0;
    static int pt_max_gain = 0;
    bool run_pt = true;
    if (Parameters::p->ls_pt_gate > 0 and pt_calls >= 2) {
        const double margin = Parameters::p->ls_pt_gate * pt_max_gain + 1.0;
        run_pt = best_solution.score_wvcp() - margin < Solution::best_score_wvcp;
    }
    if (run_pt and Parameters::p->ls_ptabu > 0 and best_solution.penalty() == 0) {
        const int before_pt = best_solution.score_wvcp();
        const auto pt_deadline =
            std::chrono::high_resolution_clock::now() +
            std::chrono::milliseconds(static_cast<long>(
                Parameters::p->ls_ptabu * 1000.0 *
                static_cast<double>(Parameters::p->max_time_local_search + 1)));
        profile_tabu(best_solution, pt_deadline);
        ++pt_calls;
        pt_max_gain = std::max(pt_max_gain, before_pt - best_solution.score_wvcp());
    }
    if (verbose) {
        print_result_ls(best_time, best_solution, turn);
    }
}

bool improve_conflicts_and_score(ProxiSolutionRedLS &solution,
                                 const int best_local_score,
                                 std::vector<bool> &tabu_list) {
    std::vector<Coloration> best_colorations;
    const int delta_wvcp = std::abs(best_local_score - solution.score_wvcp());
    int best_conflicts = 0;
    static std::vector<int> conflicting;
    solution.conflicting_vertices(conflicting);
    const bool use_atabu = atabu_on();
    const int penalty = solution.penalty();
    for (const int vertex : conflicting) {
        // for each vertices in conflict
        if (tabu_list[vertex]) {
            continue;
        }
        const int vertex_color = solution.color(vertex);
        const int *row = solution.conflicts_row(vertex);
        const int current = row[vertex_color];
        // delta_wvcp_score(vertex, color) = leave + max(0, weight - max_weight(color))
        // for every non empty colour
        const int leave = solution.delta_wvcp_score_old_color(vertex);
        const int weight = Graph::g->weights[vertex];
        for (const auto &color : solution.non_empty_colors()) {
            // for each color different from the color of the vertex
            if (color == vertex_color) {
                continue;
            }
            const int delta_conflicts = row[color] - current;
            if (delta_conflicts >= 0 or delta_conflicts > best_conflicts) {
                continue;
            }
            const int max_w = solution.max_weight(color);
            if (leave + (weight > max_w ? weight - max_w : 0) >= delta_wvcp) {
                continue;
            }
            if (use_atabu and penalty + delta_conflicts != 0 and atabu(vertex, color)) {
                continue;
            }
            // pick the move that reduce the most the number of conflicts while
            // improving the score
            if (delta_conflicts < best_conflicts) {
                best_conflicts = delta_conflicts;
                best_colorations.clear();
            }
            best_colorations.emplace_back(Coloration{vertex, color});
        }
    }

    if (best_colorations.empty()) {
        return false;
    }
    // pick one of the best improving move randomly and apply it
    const auto best_move = rd::choice(best_colorations);
    if (use_atabu) {
        atabu_set(best_move.vertex, solution.color(best_move.vertex),
                  static_cast<int>(conflicting.size()));
    }
    solution.delete_from_color(best_move.vertex);
    solution.add_to_color(best_move.vertex, best_move.color);
    // the vertex become tabu while its neighbors get out of the tabu list
    tabu_list[best_move.vertex] = true;
    for (const auto &neighbor : Graph::g->neighborhood[best_move.vertex]) {
        tabu_list[neighbor] = false;
    }

    return true;
}

bool improve_conflicts(ProxiSolutionRedLS &solution,
                       const bool with_conf,
                       std::vector<bool> &tabu_list) {
    // find a move that improve the number of conflict and possibly the score
    // if with_conf is true then the vertex must not be in the tabu list
    // if with_conf is true then the chosen vertex will be tabu after the move
    // only the unique heaviest vertex of a colour can decrease the score when it moves
    // (delta_wvcp_score < 0 requires it): scan these candidates only, in increasing
    // vertex order (same candidate list as a scan over all vertices)
    std::vector<Coloration> best_colorations;
    static std::vector<int> candidates;
    candidates.clear();
    for (const auto &color : solution.non_empty_colors()) {
        const auto &vertices = solution.colors_vertices(color);
        auto it = vertices.begin();
        const int heaviest = *it;
        ++it;
        if (it == vertices.end() or
            Graph::g->weights[*it] < Graph::g->weights[heaviest]) {
            candidates.push_back(heaviest);
        }
    }
    std::sort(candidates.begin(), candidates.end());
    for (const int vertex : candidates) {
        if (with_conf and tabu_list[vertex]) {
            continue;
        }
        for (const auto &color : solution.non_empty_colors()) {
            if (color == solution.color(vertex) or
                solution.delta_conflicts(vertex, color) > 0 or
                solution.delta_wvcp_score(vertex, color) >= 0) {
                continue;
            }
            best_colorations.emplace_back(Coloration{vertex, color});
        }
    }

    if (best_colorations.empty()) {
        return false;
    }

    const auto chosen_one = rd::choice(best_colorations);
    solution.delete_from_color(chosen_one.vertex);
    solution.add_to_color(chosen_one.vertex, chosen_one.color);

    if (with_conf) {
        tabu_list[chosen_one.vertex] = true;
    }

    return true;
}

bool ejection_move(ProxiSolutionRedLS &solution,
                   const int best_local_score,
                   std::vector<bool> &tabu_list) {
    struct Ejection {
        int v, b, u, c;
    };
    const int slack = best_local_score - solution.score_wvcp();
    const auto &weights = Graph::g->weights;
    static std::vector<int> conflicting;
    static std::vector<int> count;
    static std::vector<int> single;
    static std::vector<int> touched;
    if (count.empty()) {
        count.assign(static_cast<size_t>(Graph::g->nb_vertices) + 1, 0);
        single.assign(static_cast<size_t>(Graph::g->nb_vertices) + 1, -1);
    }
    std::vector<Ejection> best_moves;
    int best_gain = 0; // most negative change of the weighted conflicts
    solution.conflicting_vertices(conflicting);
    for (const int v : conflicting) {
        if (tabu_list[v]) {
            continue;
        }
        const int a = solution.color(v);
        const int current = solution.conflicts_row(v)[a];
        const int leave_v = solution.delta_wvcp_score_old_color(v);
        // classes b != a where v has exactly one neighbour u: from the neighbourhood when it
        // is shorter than the list of classes, else from the (count, id sum) table
        touched.clear();
        const auto &nbh = Graph::g->neighborhood[v];
        const auto &classes = solution.non_empty_colors();
        if (nbh.size() <= classes.size()) {
            for (const int neighbor : nbh) {
                const int c = solution.color(neighbor);
                if (c == a) {
                    continue;
                }
                if (count[c]++ == 0) {
                    touched.push_back(c);
                    single[c] = neighbor;
                }
            }
            size_t keep = 0;
            for (const int b : touched) {
                if (count[b] == 1) {
                    touched[keep++] = b;
                }
                count[b] = 0;
            }
            touched.resize(keep);
        } else {
            const int *nr = solution.nbr_row(v);
            for (const int b : classes) {
                if (b != a and nr[2 * b] == 1) {
                    touched.push_back(b);
                    single[b] = nr[2 * b + 1];
                }
            }
        }
        for (const int b : touched) {
            const int u = single[b];
            if (tabu_list[u]) {
                continue;
            }
            const int *row_u = solution.conflicts_row(u);
            const int gain = -current - row_u[b];
            if (gain > best_gain) {
                continue;
            }
            const int max_b = solution.max_weight(b);
            // max weight of b without u, then with v
            const int max_b_wo_u = max_b + solution.delta_wvcp_score_old_color(u);
            const int new_max_b = std::max(max_b_wo_u, weights[v]);
            const int partial = leave_v + new_max_b - max_b;
            if (partial >= slack) {
                continue;
            }
            for (const int c : solution.non_empty_colors()) {
                if (c == a or c == b or row_u[c] != 0) {
                    continue;
                }
                const int max_c = solution.max_weight(c);
                const int delta = partial + (weights[u] > max_c ? weights[u] - max_c : 0);
                if (delta >= slack) {
                    continue;
                }
                if (gain < best_gain) {
                    best_gain = gain;
                    best_moves.clear();
                }
                best_moves.push_back(Ejection{v, b, u, c});
            }
        }
    }
    if (best_moves.empty()) {
        return false;
    }
    const Ejection move = rd::choice(best_moves);
    solution.delete_from_color(move.u);
    solution.add_to_color(move.u, move.c);
    solution.delete_from_color(move.v);
    solution.add_to_color(move.v, move.b);
    for (const int neighbor : Graph::g->neighborhood[move.v]) {
        tabu_list[neighbor] = false;
    }
    for (const int neighbor : Graph::g->neighborhood[move.u]) {
        tabu_list[neighbor] = false;
    }
    tabu_list[move.v] = true;
    tabu_list[move.u] = true;
    return true;
}

bool solve_one_conflict_preserve_score(ProxiSolutionRedLS &solution,
                                       const int best_local_score,
                                       std::vector<bool> &tabu_list) {
    // select all non tabu vertices in conflict that keep the wvcp
    // score lower than the best found score
    // the vertex will be tabu after the move
    const int delta_wvcp = best_local_score - solution.score_wvcp();
    std::vector<int> vertices;
    static std::vector<int> conflicting;
    solution.conflicting_vertices(conflicting);
    for (const int vertex : conflicting) {
        if (tabu_list[vertex] == false and
            solution.delta_wvcp_score(vertex, -1) < delta_wvcp) {
            vertices.emplace_back(vertex);
        }
    }

    if (vertices.empty()) {
        return false;
    }

    const int vertex = rd::choice(vertices);
    if (atabu_on()) {
        atabu_set(vertex, solution.color(vertex), static_cast<int>(conflicting.size()));
    }
    solution.delete_from_color(vertex);
    solution.add_to_color(vertex, -1);
    tabu_list[vertex] = true;

    return true;
}

void move_heaviest_vertices(ProxiSolutionRedLS &solution) {
    // look for the best set of heaviest vertices in each colors
    // to relocate them in an other color without increasing too
    // much the number of conflicts
    float best_ratio = 0;
    int best_color = -1;
    std::set<int> best_heaviest_vertices;
    for (const auto &color1 : solution.non_empty_colors()) {
        const int max_weight1 = solution.max_weight(color1);
        int second_max = 0;
        std::set<int> heaviest_vertices;
        for (const auto &vertex : solution.colors_vertices(color1)) {
            const int vertex_weight = Graph::g->weights[vertex];
            if (vertex_weight == max_weight1) {
                heaviest_vertices.insert(vertex);
            } else if (vertex_weight > second_max) {
                second_max = vertex_weight;
            }
        }
        const int delta_move = max_weight1 - second_max;
        for (const auto &color2 : solution.non_empty_colors()) {
            if (color1 == color2) {
                continue;
            }
            int delta_conflict = 0;
            for (const auto &vertex : heaviest_vertices) {
                delta_conflict += solution.conflicts_colors(color2, vertex);
            }
            const int max_weight2 = solution.max_weight(color2);
            const int delta_wvcp = max_weight1 > max_weight2
                                       ? delta_move - max_weight1 + max_weight2
                                       : delta_move;

            float ratio = delta_conflict != 0 ? static_cast<float>(delta_wvcp) /
                                                    static_cast<float>(delta_conflict)
                                              : static_cast<float>(delta_wvcp) * 2;
            if (ratio > best_ratio or best_heaviest_vertices.empty()) {
                best_ratio = ratio;
                best_heaviest_vertices = heaviest_vertices;
                best_color = color2;
            }
        }
    }

    for (const auto &vertex : best_heaviest_vertices) {
        solution.delete_from_color(vertex);
        solution.add_to_color(vertex, best_color);
        // Probably not tabu
        // tabu_list[vertex] = 1;
    }
}

void solve_one_conflict(ProxiSolutionRedLS &solution,
                        const int best_local_score,
                        std::vector<bool> &tabu_list) {
    // Pick a random pair of vertices in conflicts and find the best move
    // possibly improving the score and the number of vertices in conflicts
    // otherwise a random move in a random color
    const int delta_wvcp = best_local_score - solution.score_wvcp();
    Coloration best_coloration{-1, -1};
    const auto &[v1, v2] = rd::choice(solution.conflict_edges());
    int best_score_conflicts = 0;
    const bool use_atabu = atabu_on();
    const int penalty = solution.penalty();
    for (const int &vertex : {v1, v2}) {
        for (const auto &color : solution.non_empty_colors()) {
            if (color == solution.color(vertex) or
                solution.delta_wvcp_score(vertex, color) >= delta_wvcp) {
                continue;
            }
            int delta_conflicts_score = solution.delta_conflicts(vertex, color);
            if (use_atabu and penalty + delta_conflicts_score != 0 and
                atabu(vertex, color)) {
                continue;
            }
            if (delta_conflicts_score < best_score_conflicts or
                best_coloration.vertex == -1) {
                best_score_conflicts = delta_conflicts_score;
                best_coloration = Coloration{vertex, color};
            }
        }
    }

    if (best_coloration.vertex == -1) {
        const std::vector<int> edges = {v1, v2};
        const int vertex = rd::choice(edges);
        std::vector<int> possible_colors{-1};
        for (const auto &color : solution.non_empty_colors()) {
            if (color != solution.color(vertex)) {
                possible_colors.emplace_back(color);
            }
        }
        const int color = rd::choice(possible_colors);

        best_coloration = Coloration{vertex, color};
    }

    if (use_atabu) {
        static std::vector<int> conflicting;
        solution.conflicting_vertices(conflicting);
        atabu_set(best_coloration.vertex, solution.color(best_coloration.vertex),
                  static_cast<int>(conflicting.size()));
    }
    solution.delete_from_color(best_coloration.vertex);
    solution.add_to_color(best_coloration.vertex, best_coloration.color);
    tabu_list[best_coloration.vertex] = true;
}
