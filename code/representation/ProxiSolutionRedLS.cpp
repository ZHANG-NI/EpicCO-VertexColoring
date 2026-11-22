#include "ProxiSolutionRedLS.hpp"

#include <algorithm>
#include <cassert>

ProxiSolutionRedLS::ProxiSolutionRedLS(Solution solution)
    : _solution(solution),
      _cap(Graph::g->nb_vertices + 1),
      _conflicts_colors(static_cast<size_t>(Graph::g->nb_vertices) *
                            (Graph::g->nb_vertices + 1),
                        0),
      _edge_weights(static_cast<size_t>(Graph::g->nb_vertices) * Graph::g->nb_vertices,
                    0),
      _nbr(2 * static_cast<size_t>(Graph::g->nb_vertices) * (Graph::g->nb_vertices + 1), 0) {
    const int n = Graph::g->nb_vertices;
    const auto &cc = solution.conflicts_colors();
    for (size_t color = 0; color < cc.size(); ++color) {
        for (int vertex = 0; vertex < n; ++vertex) {
            _conflicts_colors[static_cast<size_t>(vertex) * _cap + color] =
                cc[color][vertex];
        }
    }
    // init _edge_weights
    for (const auto &[v1, v2] : Graph::g->edges_list) {
        _edge_weights[static_cast<size_t>(v1) * n + v2] = 1;
        _edge_weights[static_cast<size_t>(v2) * n + v1] = 1;
    }
    for (int vertex = 0; vertex < n; ++vertex) {
        const int color = _solution.color(vertex);
        if (color < 0) {
            continue;
        }
        for (const int neighbor : Graph::g->neighborhood[vertex]) {
            int *p = &_nbr[2 * (static_cast<size_t>(neighbor) * _cap + color)];
            ++p[0];
            p[1] += vertex;
        }
    }
}

int ProxiSolutionRedLS::add_to_color(const int vertex, const int proposed_color) {
    const int color = _solution.add_to_color(vertex, proposed_color);
    const int n = Graph::g->nb_vertices;
    int *row = &_conflicts_colors[static_cast<size_t>(vertex) * _cap];

    if (row[color] > 0) {
        // Update penalty
        _penalty += row[color];
        for (const auto &neighbor : Graph::g->neighborhood[vertex]) {
            if (_solution.color(neighbor) == color) {
                const int lower = std::min(neighbor, vertex);
                const int higher = std::max(neighbor, vertex);
                _conflict_edges.emplace_back(std::make_tuple(lower, higher));
            }
        }
    }
    // update conflicts for neighbors
    const int *ew = &_edge_weights[static_cast<size_t>(vertex) * n];
    for (const auto &neighbor : Graph::g->neighborhood[vertex]) {
        _conflicts_colors[static_cast<size_t>(neighbor) * _cap + color] += ew[neighbor];
        int *p = &_nbr[2 * (static_cast<size_t>(neighbor) * _cap + color)];
        ++p[0];
        p[1] += vertex;
    }

    return color;
}

int ProxiSolutionRedLS::delete_from_color(const int vertex) {
    const int color = _solution.color(vertex);
    const int n = Graph::g->nb_vertices;
    const int *row = &_conflicts_colors[static_cast<size_t>(vertex) * _cap];

    // Update nb of conflicts and list of conflicting edges
    if (row[color] > 0) {
        _penalty -= row[color];
        for (const auto &neighbor : Graph::g->neighborhood[vertex]) {
            if (_solution.color(neighbor) == color) {
                const int lower = std::min(neighbor, vertex);
                const int higher = std::max(neighbor, vertex);
                int index = 0;
                for (const auto &[edge1, edge2] : _conflict_edges) {
                    if (edge1 == lower and edge2 == higher) {
                        break;
                    }
                    index++;
                }
                _conflict_edges.erase(_conflict_edges.begin() + index);
            }
        }
    }

    // update conflicts for neighbors
    const int *ew = &_edge_weights[static_cast<size_t>(vertex) * n];
    for (const int neighbor : Graph::g->neighborhood[vertex]) {
        _conflicts_colors[static_cast<size_t>(neighbor) * _cap + color] -= ew[neighbor];
        int *p = &_nbr[2 * (static_cast<size_t>(neighbor) * _cap + color)];
        --p[0];
        p[1] -= vertex;
    }

    return _solution.delete_from_color(vertex);
}

void ProxiSolutionRedLS::conflicting_vertices(std::vector<int> &out) const {
    out.clear();
    for (const auto &[v1, v2] : _conflict_edges) {
        out.push_back(v1);
        out.push_back(v2);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

void ProxiSolutionRedLS::increment_edge_weights() {
    const int n = Graph::g->nb_vertices;
    for (const auto &[edge1, edge2] : _conflict_edges) {
        ++_edge_weights[static_cast<size_t>(edge1) * n + edge2];
        ++_edge_weights[static_cast<size_t>(edge2) * n + edge1];
        ++_conflicts_colors[static_cast<size_t>(edge2) * _cap + _solution.color(edge1)];
        ++_conflicts_colors[static_cast<size_t>(edge1) * _cap + _solution.color(edge2)];
    }
    _penalty += static_cast<int>(_conflict_edges.size());
}

bool ProxiSolutionRedLS::check_solution() const {
    for (int vertex = 0; vertex < Graph::g->nb_vertices; ++vertex) {
        const int color = _solution.color(vertex);
        for (const auto &neighbor : Graph::g->neighborhood[vertex]) {
            if (_solution.color(neighbor) == color) {
                const int lower = std::min(neighbor, vertex);
                const int higher = std::max(neighbor, vertex);
                int index = 0;
                for (const auto &[edge1, edge2] : _conflict_edges) {
                    if (edge1 == lower and edge2 == higher) {
                        break;
                    }
                    index++;
                }
                assert(index < static_cast<int>(_conflict_edges.size()));
            }
        }
    }
    return _solution.check_solution();
}

[[nodiscard]] int ProxiSolutionRedLS::conflicts_colors(const int &color,
                                                       const int &vertex) const {
    return _conflicts_colors[static_cast<size_t>(vertex) * _cap + color];
}

[[nodiscard]] int ProxiSolutionRedLS::penalty() const {
    return _penalty;
}

[[nodiscard]] const std::vector<std::tuple<int, int>> &
ProxiSolutionRedLS::conflict_edges() const {
    return _conflict_edges;
}

[[nodiscard]] Solution ProxiSolutionRedLS::solution() const {
    return _solution;
}

int ProxiSolutionRedLS::delta_wvcp_score(const int vertex, const int color) const {
    return _solution.delta_wvcp_score(vertex, color);
}

[[nodiscard]] const std::vector<int> &ProxiSolutionRedLS::non_empty_colors() const {
    return _solution.non_empty_colors();
}

[[nodiscard]] int ProxiSolutionRedLS::color(const int &vertex) const {
    return _solution.color(vertex);
}

[[nodiscard]] int ProxiSolutionRedLS::score_wvcp() const {
    return _solution.score_wvcp();
}

[[nodiscard]] bool ProxiSolutionRedLS::has_conflicts(const int vertex) const {
    return _solution.has_conflicts(vertex);
}

[[nodiscard]] const std::set<int> &
ProxiSolutionRedLS::colors_vertices(const int &color) const {
    return _solution.colors_vertices(color);
}

[[nodiscard]] int ProxiSolutionRedLS::max_weight(const int &color) const {
    return _solution.max_weight(color);
}

[[nodiscard]] int ProxiSolutionRedLS::first_free_vertex() const {
    return _solution.first_free_vertex();
}
