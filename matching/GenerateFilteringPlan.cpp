//
// Created by ssunah on 11/21/18.
//

#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include <queue>
#include <utility/graphoperations.h>
#include <algorithm>

void GenerateFilteringPlan::generateTSOFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                   VertexID *&order) {
    VertexID start_vertex = selectTSOFilterStartVertex(data_graph, query_graph);
    VertexID* bfs_order;
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, bfs_order);
    GraphOperations::dfsTraversal(tree, start_vertex, query_graph->getVerticesCount(), order);
    delete[] bfs_order;
}

void GenerateFilteringPlan::generateCFLFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                  VertexID *&order, int &level_count, ui *&level_offset) {
    ui query_vertices_num = query_graph->getVerticesCount();
    VertexID start_vertex = selectCFLFilterStartVertex(data_graph, query_graph);
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    level_count = -1;
    level_offset = new ui[query_vertices_num + 1];

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        if (tree[u].level_ != level_count) {
            level_count += 1;
            level_offset[level_count] = 0;
        }

        level_offset[level_count] += 1;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];

            if (tree[u].level_ == tree[u_nbr].level_) {
                if (order_index[u_nbr] < order_index[u]) {
                    tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                }
                else {
                    tree[u].fn_[tree[u].fn_count_++] = u_nbr;
                }
            }
            else if (tree[u].level_ > tree[u_nbr].level_) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
            }
            else {
                tree[u].under_level_[tree[u].under_level_count_++] = u_nbr;
            }
        }
    }

    level_count += 1;

    ui prev_value = 0;
    for (ui i = 1; i <= level_count; ++i) {
        ui temp = level_offset[i];
        level_offset[i] = level_offset[i - 1] + prev_value;
        prev_value = temp;
    }
    level_offset[0] = 0;
}

void GenerateFilteringPlan::generateDPisoFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                    VertexID *&order) {
    VertexID start_vertex = selectDPisoStartVertex(data_graph, query_graph);
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];
            if (order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
            }
            else {
                tree[u].fn_[tree[u].fn_count_++] = u_nbr;
            }
        }
    }
}

void GenerateFilteringPlan::generateCECIFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                   VertexID *&order) {
    VertexID start_vertex = selectCECIStartVertex(data_graph, query_graph);
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];
            if (u_nbr != tree[u].parent_ && order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                tree[u_nbr].fn_[tree[u_nbr].fn_count_++] = u;
            }
        }
    }
}


// l2Match
void GenerateFilteringPlan::generateL2MatchFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                   VertexID *&order) {
    VertexID start_vertex = selectL2MatchStartVertex(data_graph, query_graph);
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];
            if (u_nbr != tree[u].parent_ && order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                tree[u_nbr].fn_[tree[u_nbr].fn_count_++] = u;
            }
        }
    }
}



VertexID GenerateFilteringPlan::selectTSOFilterStartVertex(const Graph *data_graph, const Graph *query_graph) {
    auto rank_compare = [](std::pair<VertexID, double> l, std::pair<VertexID, double> r) {
        return l.second < r.second;
    };
    // Maximum priority queue.
    std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(rank_compare)> rank_queue(rank_compare);

    // Compute the ranking.
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;
        LabelID label = query_graph->getVertexLabel(query_vertex);
        ui degree = query_graph->getVertexDegree(query_vertex);
        ui frequency = data_graph->getLabelsFrequency(label);
        double rank = frequency / (double)degree;
        rank_queue.push(std::make_pair(query_vertex, rank));
    }

    // Keep the top-3.
    while (rank_queue.size() > 3) {
        rank_queue.pop();
    }

    // Pick the one with the smallest number of candidates.
    VertexID start_vertex = 0;
    ui min_candidates_num = data_graph->getGraphMaxLabelFrequency() + 1;
    while (!rank_queue.empty()) {
        VertexID query_vertex = rank_queue.top().first;

        if (rank_queue.size() == 1) {
            ui count;
            FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
            if (count < min_candidates_num) {
                start_vertex = query_vertex;
            }
        }
        else {
            LabelID label = query_graph->getVertexLabel(query_vertex);
            ui frequency = data_graph->getLabelsFrequency(label);
            if (frequency / (double)data_graph->getVerticesCount() <= 0.05) {
                ui count;
                FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
                if (count < min_candidates_num) {
                    start_vertex = query_vertex;
                    min_candidates_num = count;
                }
            }
        }
        rank_queue.pop();
    }

    return start_vertex;
}

VertexID GenerateFilteringPlan::selectCFLFilterStartVertex(const Graph *data_graph, const Graph *query_graph) {
    auto rank_compare = [](std::pair<VertexID, double> l, std::pair<VertexID, double> r) {
        return l.second < r.second;
    };

    std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(rank_compare)> rank_queue(rank_compare);

    // Compute the ranking.
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;

        if (query_graph->get2CoreSize() == 0 || query_graph->getCoreValue(query_vertex) > 1) {
            LabelID label = query_graph->getVertexLabel(query_vertex);
            ui degree = query_graph->getVertexDegree(query_vertex);
            ui frequency = data_graph->getLabelsFrequency(label);
            double rank = frequency / (double) degree;
            rank_queue.push(std::make_pair(query_vertex, rank));
        }
    }

    // Keep the top-3.
    while (rank_queue.size() > 3) {
        rank_queue.pop();
    }

    VertexID start_vertex = 0;
    double min_score = data_graph->getGraphMaxLabelFrequency() + 1;

    while (!rank_queue.empty()) {
        VertexID query_vertex = rank_queue.top().first;
        ui count;
        FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
        double cur_score = count / (double) query_graph->getVertexDegree(query_vertex);

        if (cur_score < min_score) {
            start_vertex = query_vertex;
            min_score = cur_score;
        }
        rank_queue.pop();
    }

    return start_vertex;
}

VertexID GenerateFilteringPlan::selectDPisoStartVertex(const Graph *data_graph, const Graph *query_graph) {
    double min_score = data_graph->getVerticesCount();
    VertexID start_vertex = 0;

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        ui degree = query_graph->getVertexDegree(i);
        if (degree <= 1)
            continue;

        ui count = 0;
        FilterVertices::computeCandidateWithLDF(data_graph, query_graph, i, count);
        double cur_score = count / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }

    return start_vertex;
}

VertexID GenerateFilteringPlan::selectCECIStartVertex(const Graph *data_graph, const Graph *query_graph) {
    double min_score = data_graph->getVerticesCount();
    VertexID start_vertex = 0;

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        ui degree = query_graph->getVertexDegree(i);
        ui count = 0;
        FilterVertices::computeCandidateWithNLF(data_graph, query_graph, i, count);
        double cur_score = count / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }

    return start_vertex;
}


// l2Match
VertexID GenerateFilteringPlan::selectL2MatchStartVertex(const Graph *data_graph, const Graph *query_graph) {
    double min_score = data_graph->getVerticesCount();
    VertexID start_vertex = 0;

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        ui degree = query_graph->getVertexDegree(i);
        ui count = 0;
        FilterVertices::computeCandidateWithNLF_LPF(data_graph, query_graph, i, count);
        double cur_score = count / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }

    return start_vertex;
}

// l2Match with Spatial Ordering
void GenerateFilteringPlan::generateL2MatchFilterPlanB(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                  VertexID *&order) {
    ui query_vertices_num = query_graph->getVerticesCount();
    VertexID start_vertex = selectL2MatchStartVertex(data_graph, query_graph);
    ui* spatial_index_ = new ui[query_vertices_num];
    std::vector<bool> visited_query_vertex(query_vertices_num, false);
    std::queue<ui> q;
    
    q.push(start_vertex);
    visited_query_vertex[start_vertex] = true;

    tree = new TreeNode[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        tree[i].initialize(query_vertices_num);
    }
    order = new VertexID[query_vertices_num];
    tree[start_vertex].level_ = 0;
    tree[start_vertex].id_ = start_vertex;

    ui si = 0;
    spatial_index_[start_vertex] = si;
    ui last_v = start_vertex;
    std::vector<std::vector<ui>> similar_si_nbrs(query_vertices_num); // neighbors of similar spatial index
    std::vector<std::vector<ui>> diff_si_nbrs(query_vertices_num); // neighbors of different spatial index
    
    // BFS Traversal to calculate spatial difference
    while (!q.empty()) {
        VertexID u = q.front();
        q.pop();
        // all pair of connected vertices either has:
        //      difference of -1 (previous depth) or +1 (next depth),
        //      or same spatial index (same level of depth)
        
        ui deg;
        const ui* adjs = query_graph->getVertexNeighbors(u, deg);
        
        for (ui j=0; j<deg; ++j) {
            ui v = adjs[j];
            if (!visited_query_vertex[v]) {
                visited_query_vertex[v] = true;
                q.push(v);
                tree[v].id_ = v;
                tree[v].parent_ = u;
                tree[v].level_ = tree[u].level_ + 1;
                tree[u].children_[tree[u].children_count_++] = v;
                spatial_index_[v] = spatial_index_[u] + 1;
            }
            else if (spatial_index_[v] == spatial_index_[u]) {
                similar_si_nbrs[v].emplace_back(u);
            } else {
                diff_si_nbrs[v].emplace_back(u);
            }
        }

        if (u == last_v) {
            si += 1;
            // let the last vertex in next depth = last element in the queue
            last_v = q.back();
        }
    }

    for (ui i=0; i<query_vertices_num; ++i) {
        ui si = spatial_index_[i];
        // O(|E|)
        // sort neighbors of same spatial index
        std::sort(similar_si_nbrs[i].begin(), similar_si_nbrs[i].end(), [&similar_si_nbrs](const ui& l, const ui& r){
            return similar_si_nbrs[l].size() < similar_si_nbrs[r].size();
        });
        // sort neighbors of different spatial index
        std::sort(diff_si_nbrs[i].begin(), diff_si_nbrs[i].end(), [&similar_si_nbrs](const ui& l, const ui& r){
            return similar_si_nbrs[l].size() < similar_si_nbrs[r].size();
        });
    }

    // BFS (NTE) Traversal to generate order
    ui count = 0;
    q.push(start_vertex);
    std::fill(visited_query_vertex.begin(), visited_query_vertex.end(), false);
    visited_query_vertex[start_vertex] = true;
    std::vector<ui> deg_one_nbrs;

    while (!q.empty()) {
        VertexID u = q.front();
        q.pop();
        order[count++] = u;

        ui deg;
        ui nbrs_count = 0;
        const ui* adjs = query_graph->getVertexNeighbors(u, deg);
        
        // prioritize neighbor with similar spatial index (same breath)
        while (!similar_si_nbrs[u].empty()) {
            nbrs_count += 1;

            // prioritize neighbor with maximum neighbors of in same depth
            ui v = similar_si_nbrs[u].back();
            similar_si_nbrs[u].pop_back();
            
            if (!visited_query_vertex[v]) {
                visited_query_vertex[v] = true;
                
                // if a neighbor is at same depth, then this neighbor is connected to u and a vertex in previous depth
                //      thus, the degree of this neighbor must be >= 2. No checking for leaf vertex is required
                //if (query_graph->getVertexDegree(v) == 1) {
                //    deg_one_nbrs.push_back(v);
                //    continue;
                //}
                q.push(v);
            }
        }

        // prioritize forward neighbor in next depth with most connection within its depth
        while (!diff_si_nbrs[u].empty()) {
            ui v = diff_si_nbrs[u].back();
            diff_si_nbrs[u].pop_back();

            if (!visited_query_vertex[v]) {
                visited_query_vertex[v] = true;
                if (query_graph->getVertexDegree(v) == 1) {
                    deg_one_nbrs.push_back(v);
                    continue;
                }
                q.push(v);
            }
        }
    }
    
    // finally, append degree one (leaf) vertex to order
    for (ui u : deg_one_nbrs) {
        order[count++] = u;
    }

    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    // compute backward, forward neighbors and parent
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];
            if (u_nbr != tree[u].parent_ && order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                tree[u_nbr].fn_[tree[u_nbr].fn_count_++] = u;
            }
        }
    }
}