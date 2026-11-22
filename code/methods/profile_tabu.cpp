#include "profile_tabu.hpp"

#include <algorithm>
#include <limits>
#include <numeric>

#include "../utils/random_generator.hpp"

namespace {

bool time_over(const std::chrono::high_resolution_clock::time_point &deadline) {
    return std::chrono::high_resolution_clock::now() >= deadline or
           Parameters::p->time_limit_reached();
}

/** legal descent: move a vertex to another class without conflict when it lowers the score */
void legal_descent(Solution &solution) {
    const int n = Graph::g->nb_vertices;
    bool improved = true;
    while (improved) {
        improved = false;
        for (int v = 0; v < n; ++v) {
            const int cv = solution.color(v);
            if (solution.delta_wvcp_score_old_color(v) >= 0) {
                continue; // leaving its class does not lower the score
            }
            int best_c = -2;
            int best_d = 0;
            for (const int c : solution.non_empty_colors()) {
                if (c == cv or solution.conflicts_colors(c, v) != 0) {
                    continue;
                }
                const int d = solution.delta_wvcp_score(v, c);
                if (d < best_d) {
                    best_d = d;
                    best_c = c;
                }
            }
            if (best_c != -2) {
                solution.delete_from_color(v);
                solution.add_to_color(v, best_c);
                improved = true;
            }
        }
    }
}

struct TabuCol {
    int n;
    int k;
    std::vector<int> level;
    std::vector<int> col;
    std::vector<int> gamma; // k*n, class-major: gamma[c * n + v]
    std::vector<long> tabu; // n*k
    std::vector<int> conf_list;
    std::vector<int> conf_pos;
    int nb_conf = 0;

    void add_conf(int v) {
        conf_pos[v] = static_cast<int>(conf_list.size());
        conf_list.push_back(v);
    }
    void del_conf(int v) {
        const int p = conf_pos[v];
        const int last = conf_list.back();
        conf_list[p] = last;
        conf_pos[last] = p;
        conf_list.pop_back();
        conf_pos[v] = -1;
    }

    void init() {
        gamma.assign(static_cast<size_t>(n) * k, 0);
        tabu.assign(static_cast<size_t>(n) * k, 0);
        conf_pos.assign(n, -1);
        conf_list.clear();
        nb_conf = 0;
        for (const auto &[u, v] : Graph::g->edges_list) {
            ++gamma[static_cast<size_t>(col[v]) * n + u];
            ++gamma[static_cast<size_t>(col[u]) * n + v];
            if (col[u] == col[v])
                ++nb_conf;
        }
        for (int v = 0; v < n; ++v) {
            if (gamma[static_cast<size_t>(col[v]) * n + v] > 0)
                add_conf(v);
        }
    }

    void move(int v, int c) {
        const int old = col[v];
        int *go = &gamma[static_cast<size_t>(old) * n];
        int *gc = &gamma[static_cast<size_t>(c) * n];
        nb_conf += gc[v] - go[v];
        col[v] = c;
        for (const int u : Graph::g->neighborhood[v]) {
            --go[u];
            ++gc[u];
            if (col[u] == old) {
                if (go[u] == 0)
                    del_conf(u);
            } else if (col[u] == c) {
                if (gc[u] == 1)
                    add_conf(u);
            }
        }
        const bool was = conf_pos[v] >= 0;
        const bool is = gc[v] > 0;
        if (was and not is)
            del_conf(v);
        else if (not was and is)
            add_conf(v);
    }

    /** returns true when a conflict-free assignment was found */
    bool run(const long stall,
             const std::chrono::high_resolution_clock::time_point &deadline) {
        int best_conf = nb_conf;
        long iter = 0;
        long iter_best = 0;
        std::uniform_int_distribution<int> tdis(0, 9);
        // allowed classes (level >= w(v)) per threshold, in increasing class index: the scan
        // below visits exactly the classes the full scan would accept, in the same order
        std::vector<int> values(level.begin(), level.end());
        std::sort(values.begin(), values.end(), std::greater<int>());
        values.erase(std::unique(values.begin(), values.end()), values.end());
        const int d = static_cast<int>(values.size());
        std::vector<int> start(static_cast<size_t>(d) + 2, 0);
        std::vector<int> flat;
        flat.reserve(static_cast<size_t>(d) * k);
        start[0] = 0; // list 0 = empty (no class allowed)
        for (int t = 0; t < d; ++t) {
            for (int c = 0; c < k; ++c) {
                if (level[c] >= values[t])
                    flat.push_back(c);
            }
            start[t + 2] = static_cast<int>(flat.size());
        }
        start[1] = 0;
        // list of vertex v: [start[thr[v]], start[thr[v] + 1]) with thr = 1 + last t, values[t] >= w
        std::vector<int> thr(n, 0);
        for (int v = 0; v < n; ++v) {
            const int wv = Graph::g->weights[v];
            // number of distinct values >= wv
            const int cnt = static_cast<int>(
                std::upper_bound(values.begin(), values.end(), wv, std::greater<int>()) -
                values.begin());
            thr[v] = cnt; // cnt = 0 -> empty list [start[0], start[1])
        }
        while (nb_conf > 0) {
            ++iter;
            if (iter - iter_best > stall)
                return false;
            if ((iter & 255) == 0 and time_over(deadline))
                return false;
            int best_delta = std::numeric_limits<int>::max();
            int bv = -1;
            int bc = -1;
            int nb_ties = 0;
            for (const int v : conf_list) {
                const int *g = &gamma[v];
                const long *t = &tabu[static_cast<size_t>(v) * k];
                const int cv = col[v];
                const int cur = g[static_cast<size_t>(cv) * n];
                const int *lb = flat.data() + start[thr[v]];
                const int *le = flat.data() + start[thr[v] + 1];
                for (const int *it = lb; it != le; ++it) {
                    const int c = *it;
                    if (c == cv)
                        continue;
                    const int delta = g[static_cast<size_t>(c) * n] - cur;
                    if (delta > best_delta)
                        continue;
                    if (t[c] > iter and nb_conf + delta >= best_conf)
                        continue;
                    if (delta < best_delta) {
                        best_delta = delta;
                        bv = v;
                        bc = c;
                        nb_ties = 1;
                    } else {
                        ++nb_ties;
                        if (std::uniform_int_distribution<int>(0, nb_ties - 1)(rd::generator) == 0) {
                            bv = v;
                            bc = c;
                        }
                    }
                }
            }
            if (bv == -1) {
                continue; // every move tabu
            }
            const int old = col[bv];
            move(bv, bc);
            tabu[static_cast<size_t>(bv) * k + old] =
                iter + tdis(rd::generator) +
                static_cast<long>(0.6 * static_cast<double>(conf_list.size()));
            if (nb_conf < best_conf) {
                best_conf = nb_conf;
                iter_best = iter;
            }
        }
        return true;
    }
};

} // namespace

void profile_tabu(Solution &solution,
                  const std::chrono::high_resolution_clock::time_point &deadline) {
    const int n = Graph::g->nb_vertices;
    const auto &w = Graph::g->weights;
    const long base_stall =
        std::max(1000L, static_cast<long>(Parameters::p->ls_pt_stall * static_cast<double>(n)));
    const long max_stall = std::max(
        base_stall,
        static_cast<long>(Parameters::p->ls_pt_stall_max * static_cast<double>(n)));
    // escalating depth: every lowering is first tried with the base stall; only when all
    // of them fail and time is left, they are tried again with a 10 times larger stall
    long stall = base_stall;
    legal_descent(solution);
#ifdef PT_DEBUG
    const int s0 = solution.score_wvcp();
    const auto t0 = std::chrono::high_resolution_clock::now();
    int nb_att = 0;
#endif
    while (not time_over(deadline)) {
        // classes of the current solution, relabelled 0..k-1
        const auto &colors = solution.non_empty_colors();
        const int k = static_cast<int>(colors.size());
        std::vector<int> index(static_cast<size_t>(*std::max_element(colors.begin(), colors.end())) + 1, -1);
        for (int i = 0; i < k; ++i)
            index[colors[i]] = i;
        std::vector<int> col(n);
        std::vector<int> level(k, 0);
        std::vector<int> second(k, 0);
        for (int v = 0; v < n; ++v) {
            col[v] = index[solution.color(v)];
            level[col[v]] = std::max(level[col[v]], w[v]);
        }
        for (int v = 0; v < n; ++v) {
            if (w[v] < level[col[v]])
                second[col[v]] = std::max(second[col[v]], w[v]);
        }
        // candidate lowerings: class c from level[c] to second[c]; possible when every
        // displaced vertex has another class of level >= its weight
        std::vector<int> best_other(k, 0); // max level over the other classes
        {
            int m1 = 0, m2 = 0, arg = -1;
            for (int c = 0; c < k; ++c) {
                if (level[c] > m1) {
                    m2 = m1;
                    m1 = level[c];
                    arg = c;
                } else if (level[c] > m2) {
                    m2 = level[c];
                }
            }
            for (int c = 0; c < k; ++c)
                best_other[c] = (c == arg) ? m2 : m1;
        }
        std::vector<int> cands;
        for (int c = 0; c < k; ++c) {
            if (best_other[c] >= level[c]) // the displaced vertices (w <= level[c]) fit elsewhere
                cands.push_back(c);
        }
        // predicted difficulty of each lowering: conflicts left after placing its displaced
        // vertices greedily (exact for the initial placement: displaced vertices share a
        // class, so they are pairwise non-adjacent and do not affect each other's gamma)
        std::vector<int> est(k, 0);
        {
            std::vector<int> gam(static_cast<size_t>(n) * k, 0);
            for (const auto &[u, v] : Graph::g->edges_list) {
                ++gam[static_cast<size_t>(u) * k + col[v]];
                ++gam[static_cast<size_t>(v) * k + col[u]];
            }
            std::vector<char> is_cand(k, 0);
            for (const int c : cands)
                is_cand[c] = 1;
            for (int v = 0; v < n; ++v) {
                const int c = col[v];
                if (not is_cand[c] or w[v] <= second[c])
                    continue;
                int bg = std::numeric_limits<int>::max();
                const int *g = &gam[static_cast<size_t>(v) * k];
                for (int c2 = 0; c2 < k; ++c2) {
                    if (c2 != c and level[c2] >= w[v] and g[c2] < bg)
                        bg = g[c2];
                }
                est[c] += bg;
            }
        }
        std::shuffle(cands.begin(), cands.end(), rd::generator);
        std::stable_sort(cands.begin(), cands.end(), [&](int a, int b) {
            if (est[a] != est[b])
                return est[a] < est[b];
            return level[a] - second[a] < level[b] - second[b];
        });

        bool success = false;
        for (const int c : cands) {
            if (time_over(deadline))
                break;
            TabuCol tc;
            tc.n = n;
            tc.k = k;
            tc.level = level;
            tc.level[c] = second[c];
            tc.col = col;
            tc.init();
            // displaced vertices go to the allowed class with fewest conflicts
            for (int v = 0; v < n; ++v) {
                if (tc.col[v] != c or w[v] <= tc.level[c])
                    continue;
                int bc = -1, bg = std::numeric_limits<int>::max(), ties = 0;
                for (int c2 = 0; c2 < k; ++c2) {
                    if (c2 == c or tc.level[c2] < w[v])
                        continue;
                    const int g = tc.gamma[static_cast<size_t>(c2) * n + v];
                    if (g < bg) {
                        bg = g;
                        bc = c2;
                        ties = 1;
                    } else if (g == bg) {
                        ++ties;
                        if (std::uniform_int_distribution<int>(0, ties - 1)(rd::generator) == 0)
                            bc = c2;
                    }
                }
                tc.move(v, bc);
            }
            if (tc.run(stall, deadline)) {
                // rebuild the solution from the conflict-free assignment
                Solution s;
                std::vector<int> new_id(k, -1);
                for (int v = 0; v < n; ++v) {
                    const int cl = tc.col[v];
                    if (new_id[cl] == -1)
                        new_id[cl] = s.add_to_color(v, -1);
                    else
                        s.add_to_color(v, new_id[cl]);
                }
                while (s.first_free_vertex() < solution.first_free_vertex())
                    s.increment_first_free_vertex();
                legal_descent(s);
                if (s.penalty() == 0 and s.score_wvcp() < solution.score_wvcp()) {
                    solution = s;
                    success = true;
                }
                break;
            }
        }
        if (success) {
            stall = base_stall;
        } else if (stall < max_stall) {
            stall = std::min(max_stall, stall * 10);
        } else {
            break;
        }
    }
#ifdef PT_DEBUG
    fmt::print(stderr, "PT {} -> {} in {} ms\n", s0, solution.score_wvcp(),
               std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::high_resolution_clock::now() - t0).count());
#endif
}
