#pragma once

#include "Solution.hpp"

/**
 * @brief Proxi solution for RedLS
 *
 * Add conflicts on edges and other way to manage of conflicts on colors
 *
 */
class ProxiSolutionRedLS {

    /** @brief the current solution*/
    Solution _solution;
    /** @brief number of conflicts in the current solution*/
    int _penalty{0};

    /** @brief For each color, for each vertex, number of neighbors in the color*/
    /** @brief number of colour slots per vertex in _conflicts_colors (a colour id is < n + 1)*/
    int _cap{0};
    /** @brief vertex-major: _conflicts_colors[vertex * _cap + color] = weighted conflicts*/
    std::vector<int> _conflicts_colors{};

    /** @brief list of conflicting edges in the solution*/
    std::vector<std::tuple<int, int>> _conflict_edges{};
    /** @brief edge weights (for RedLS)*/
    /** @brief flat n x n edge weights*/
    std::vector<int> _edge_weights{};
    /** @brief vertex-major unweighted neighbourhood per colour: _nbr[2 * (vertex * _cap + color)]
     * = number of neighbours of vertex in color, [.. + 1] = sum of their ids (the neighbour itself
     * when the number is 1)*/
    std::vector<int> _nbr{};

  public:
    ProxiSolutionRedLS() = delete;
    ProxiSolutionRedLS(Solution solution);

    int add_to_color(const int vertex, const int proposed_color);

    int delete_from_color(const int vertex);

    void increment_edge_weights();

    bool check_solution() const;

    [[nodiscard]] inline int delta_conflicts(const int vertex, const int color) const {
        const int *row = &_conflicts_colors[static_cast<size_t>(vertex) * _cap];
        return row[color] - row[_solution.color(vertex)];
    }

    /** @brief pointer to the weighted conflicts of the vertex, indexed by colour*/
    [[nodiscard]] inline const int *conflicts_row(const int vertex) const {
        return &_conflicts_colors[static_cast<size_t>(vertex) * _cap];
    }

    /** @brief pointer to the (count, id sum) pairs of the neighbours of the vertex, by colour*/
    [[nodiscard]] inline const int *nbr_row(const int vertex) const {
        return &_nbr[2 * static_cast<size_t>(vertex) * _cap];
    }

    /** @brief score difference when the vertex leaves its colour*/
    [[nodiscard]] inline int delta_wvcp_score_old_color(const int vertex) const {
        return _solution.delta_wvcp_score_old_color(vertex);
    }

    /** @brief conflicting vertices in increasing order (from the conflict edges)*/
    void conflicting_vertices(std::vector<int> &out) const;

    [[nodiscard]] int conflicts_colors(const int &color, const int &vertex) const;

    [[nodiscard]] int penalty() const;

    [[nodiscard]] const std::vector<std::tuple<int, int>> &conflict_edges() const;

    [[nodiscard]] Solution solution() const;

    // Getters to original solution

    int delta_wvcp_score(const int vertex, const int color) const;

    [[nodiscard]] const std::vector<int> &non_empty_colors() const;

    [[nodiscard]] int color(const int &vertex) const;

    [[nodiscard]] int score_wvcp() const;

    [[nodiscard]] bool has_conflicts(const int vertex) const;

    [[nodiscard]] const std::set<int> &colors_vertices(const int &color) const;

    [[nodiscard]] int max_weight(const int &color) const;

    [[nodiscard]] int first_free_vertex() const;
};
