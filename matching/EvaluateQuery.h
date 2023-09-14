//
// Created by ssunah on 11/20/18.
//

#ifndef SUBGRAPHMATCHING_EVALUATEQUERY_H
#define SUBGRAPHMATCHING_EVALUATEQUERY_H

#include "graph/graph.h"
#include "utility/QFliter.h"
#include "utility/computesetintersection.h"
#include <vector>
#include <queue>
#include <bitset>
#include <cstring>
#include "utility/pretty_print.h"

// Min priority queue.
auto extendable_vertex_compare = [](std::pair<std::pair<VertexID, ui>, ui> l, std::pair<std::pair<VertexID, ui>, ui> r) {
    if (l.first.second == 1 && r.first.second != 1) {
        return true;
    }
    else if (l.first.second != 1 && r.first.second == 1) {
        return false;
    }
    else
    {
        return l.second > r.second;
    }
};

typedef std::priority_queue<std::pair<std::pair<VertexID, ui>, ui>, std::vector<std::pair<std::pair<VertexID, ui>, ui>>,
        decltype(extendable_vertex_compare)> dpiso_min_pq;

class EvaluateQuery {
public:
    static size_t exploreGraph(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                                  ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t exploreGraphJR(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                            ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t LFTJ(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                           ui *order, size_t output_limit_num, size_t &call_count);

    static size_t
    exploreGraphQLStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                            size_t output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t
    exploreQuickSIStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                            ui *pivot, size_t output_limit_num, size_t &call_count);

    static size_t exploreDPisoStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                    Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                    ui **weight_array, ui *order, size_t output_limit_num,
                                    size_t &call_count, TimeOutException* timeout);

    static size_t exploreDPisoRecursiveStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                             Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                             ui **weight_array, ui *order, size_t output_limit_num,
                                             size_t &call_count);

    static size_t exploreCECIStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                      ui *candidates_count,
                                      std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                      std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                      ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t exploreL2Match(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                      ui *candidates_count,
                                      std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                      std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                      ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t exploreL2MatchB(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                      ui *candidates_count,
                                      std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                      ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t exploreL2MatchC(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                      ui *candidates_count,
                                      std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                      std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                      ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t exploreL2MatchNew(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, TreeNode *tree, ui **candidates,
                                  ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count, TimeOutException* timeout);

    static size_t exploreUllman(const Graph* data_graph, const Graph* query_graph, size_t output_limit_num, TimeOutException* timeout);
    static bool recurseUllman(bool* visited_candidates, ui depth, const Graph* data_graph, const Graph* query_graph, ui** mappings, ui& call_count, size_t& embedding_count, size_t output_limit_num, TimeOutException* timeout);
    static void pruneUllman(const Graph* data_graph, const Graph* query_graph, ui** mappings);

#if ENABLE_QFLITER == 1
    static BSRGraph*** qfliter_bsr_graph_;
    static int* temp_bsr_base1_;
    static int* temp_bsr_state1_;
    static int* temp_bsr_base2_;
    static int* temp_bsr_state2_;
#endif

    static bool exit_;

#ifdef DISTRIBUTION
    static size_t* distribution_count_;
#endif
private:
    static void generateBN(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count);
    static void generateBN(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count);
    static void generateTree(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count,
                                int*& childs, int*& jump);
    static void generateTree(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count,
                            int* childs, int* jump);
    static void allocateBuffer(const Graph *query_graph, const Graph *data_graph, ui *candidates_count, ui *&idx,
                                   ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                                   ui **&valid_candidate_idx, bool *&visited_vertices);
    static void releaseBuffer(ui query_vertices_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn, ui *bn_count);

    static void generateValidCandidateIndex(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                            ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                            bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot,
                                            ui **candidates);

    static void generateValidCandidateIndexJR(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                            ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                            ui **bn, ui *bn_cnt, ui *order, ui *pivot, ui **candidates);

    static void generateValidCandidateIndexJRNew(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                            ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                            ui **bn, ui *bn_cnt, ui *order, ui pivot, ui **candidates);

    static void generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order, ui *&temp_buffer);

    static void generateValidCandidates(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                        bool* visited_vertices, ui **bn, ui *bn_cnt, ui* order, ui **candidates, ui* candidates_count);

    static void generateValidCandidates(const Graph *query_graph, const Graph *data_graph, ui depth, ui *embedding,
                                            ui *idx_count, ui **valid_candidate, bool *visited_vertices, ui **bn, ui *bn_cnt,
                                            ui *order, ui *pivot);
    static void generateValidCandidates(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            ui *&temp_buffer, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);

    static void generateValidCandidatesL2MatchB(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order, TreeNode *tree,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);    

    static void generateValidCandidatesL2MatchC(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            TreeNode *tree, bool *visited_vertices,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);  

    static void updateExtendableVertex(ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                          Edges ***edge_matrix, ui *&temp_buffer, ui **weight_array,
                                          TreeNode *tree, VertexID mapped_vertex, ui *extendable,
                                          std::vector<dpiso_min_pq> &vec_rank_queue, const Graph *query_graph);

    static void restoreExtendableVertex(TreeNode* tree, VertexID unmapped_vertex, ui *extendable);
    static void generateValidCandidateIndex(VertexID vertex, ui *idx_embedding, ui *idx_count, ui *&valid_candidate_index,
                                            Edges ***edge_matrix, ui *bn, ui bn_cnt, ui *&temp_buffer);

    static void computeAncestor(const Graph *query_graph, TreeNode *tree, VertexID *order,
                                std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors);

    static void computeAncestor(const Graph *query_graph, ui** bn, ui* bn_cnt, VertexID *order,
                                std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors);

    static void computeAncestor(const Graph *query_graph, VertexID *order, std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors);

    static std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> exploreDPisoBacktrack(ui max_depth, ui depth, VertexID mapped_vertex, TreeNode *tree, ui *idx_embedding,
                                                     ui *embedding, std::unordered_map<VertexID, VertexID> &reverse_embedding,
                                                     bool *visited_vertices, ui *idx_count, ui **valid_candidate_index,
                                                     Edges ***edge_matrix,
                                                     std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors,
                                                     dpiso_min_pq rank_queue, ui **weight_array, ui *&temp_buffer, ui *extendable,
                                                     ui **candidates, size_t &embedding_count, size_t &call_count,
                                                     const Graph *query_graph);
};

#include "utility/pretty_print.h"

#if ENABLE_QFLITER == 1
BSRGraph ***EvaluateQuery::qfliter_bsr_graph_;
int *EvaluateQuery::temp_bsr_base1_ = nullptr;
int *EvaluateQuery::temp_bsr_state1_ = nullptr;
int *EvaluateQuery::temp_bsr_base2_ = nullptr;
int *EvaluateQuery::temp_bsr_state2_ = nullptr;
#endif

#ifdef SPECTRUM
bool EvaluateQuery::exit_;
#endif

#ifdef DISTRIBUTION
size_t* EvaluateQuery::distribution_count_;
#endif

void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr] && nbr != pivot[i]) {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        #if DEBUGGING_MODE == 1
        std::cout << "\tcurr u =\t" << vertex << std::endl;
        #endif

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
                #if DEBUGGING_MODE == 1
                std::cout << "\t\tparent of u[" << vertex << "] IS v[" << nbr << "](" << (bn_count[i]) << ")" << std::endl;
                #endif
            }
        }

        visited_vertices[vertex] = true;
    }
}

size_t
EvaluateQuery::exploreGraph(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                            ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count, TimeOutException* timeout) {
    // Generate the bn.
    int max_depth = query_graph->getVerticesCount();
#ifdef JR
    // recompute fn/child candidate set if bn/parent's candidate has changed
    int* childs = new int[max_depth]; // fn and child of any query vertex
    int redo = 0; // flag of query vertices that require recomputation of candidate set
    int* jump = new int[max_depth]; // closest bn/parent index to jump to if enumeration failed at any query vertex
#endif
    ui **bn;
    ui *bn_count;

#ifdef JR
    generateTree(query_graph, order, pivot, bn, bn_count, childs, jump);
#else
    generateBN(query_graph, order, pivot, bn, bn_count);
#endif

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            idx[cur_depth] += 1;

#ifdef JR
            if (visited_vertices[v]) {
                // find next non-visited valid candidate
                continue;
            }
            // REDO
            redo |= childs[cur_depth];
#endif

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                cur_depth += 1;
                idx[cur_depth] = 0;
#ifdef JR
                if (redo & (1 << cur_depth)) {
                    redo -= (1 << cur_depth);

                    call_count += 1;
                    idx_count[cur_depth] = 0;
                    generateValidCandidateIndex(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx,
                                            edge_matrix, visited_vertices, bn, bn_count, order, pivot, candidates);
                }

                if (idx_count[cur_depth] == 0) {
                    int jump_depth = jump[cur_depth];
                    for (int i=cur_depth-1; i>= jump_depth; --i) {
                        ui v = embedding[order[i]];
                        visited_vertices[v] = false; // unvisit and unmap previous query vertex's candidate
                    }
                    cur_depth = jump_depth; // jump to closest bn/parent
                }
#else
                call_count += 1;
                generateValidCandidateIndex(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx,
                                            edge_matrix, visited_vertices, bn, bn_count, order, pivot, candidates);
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }


    // Release the buffer.
    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);
#ifdef JR
    delete[] childs;
    delete[] jump;
#endif
    return embedding_cnt;
}

void EvaluateQuery::generateTree(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count,
                                    int*& childs, int*& jump) {
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    std::vector<int> order_idx(query_vertices_num);

    for (ui i = 0; i < query_vertices_num; ++i) {
        bn[i] = new ui[query_vertices_num];
        order_idx[order[i]] = i;
        childs[i] = 0;
        jump[i] = 0;
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;

    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        int max_idx = 0;
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];
            if (!visited_vertices[nbr])
                continue;
            if (nbr != pivot[i]) {
                bn[i][bn_count[i]++] = nbr;

            }
            childs[ order_idx[nbr] ] += (1 << i);
            if (max_idx < order_idx[nbr])
                max_idx = order_idx[nbr];

        }

        jump[i] = max_idx;
        visited_vertices[vertex] = true;
    }
}

void EvaluateQuery::generateTree(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count, int* childs, int* jump) {
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    std::vector<int> order_idx(query_vertices_num);

    for (ui i = 0; i < query_vertices_num; ++i) {
        bn[i] = new ui[query_vertices_num];
        order_idx[order[i]] = i;
        childs[i] = 0;
        jump[i] = 0;
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        
        int max_idx = 0;
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
                
                childs[ order_idx[nbr] ] += (1 << i);
                if (max_idx < order_idx[nbr])
                    max_idx = order_idx[nbr];
            }
        }

        jump[i] = max_idx;
        visited_vertices[vertex] = true;
    }
}

// l2Match with Jump & Redo
size_t
EvaluateQuery::exploreGraphJR(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                            ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count, TimeOutException* timeout) {
    // Jump & Redo variables
    int max_depth = query_graph->getVerticesCount();
    // recompute fn/child candidate set if bn/parent's candidate has changed
    int* childs = new int[max_depth]; // fn and child of any query vertex
    int redo = 0; // flag of query vertices that require recomputation of candidate set
    int* jump = new int[max_depth]; // closest bn/parent index to jump to if enumeration failed at any query vertex
    //size_t* partial_embedding_count = new size_t[max_depth]; // embedding count since traversal at any query vertex
    // Jump & Redo variables end
    
    // Generate the bn.
    ui **bn;
    ui *bn_count;
    generateTree(query_graph, order, pivot, bn, bn_count, childs, jump);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    //int *reverse_embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

    //reverse_embedding = new int[data_graph->getVerticesCount()];


    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {
                // find next non-visited valid candidate
                continue;
            }

            redo |= childs[cur_depth];

            embedding[u] = v;
            //reverse_embedding[v] = cur_depth;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                //reverse_embedding[v] = 0;
                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                cur_depth += 1;
                idx[cur_depth] = 0;
                //ui u_f = order[cur_depth];
                //partial_embedding_count[cur_depth] = embedding_cnt;

                if (redo & (1 << cur_depth)) {
                    redo -= (1 << cur_depth);

                    call_count += 1;
                    idx_count[cur_depth] = 0;
                    generateValidCandidateIndexJR(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                                valid_candidate_idx,
                                                edge_matrix, bn, bn_count, order, pivot, candidates);
                }

                if (idx_count[cur_depth] == 0) {
                    int jump_depth = jump[cur_depth];
                    for (int i=cur_depth-1; i>= jump_depth; --i) {
                        ui v = embedding[order[i]];
                        //reverse_embedding[v] = 0;
                        visited_vertices[v] = false; // unvisit and unmap previous query vertex's candidate
                    }
                    //visited_vertices[ embedding[ order[jump_depth] ] ] = false;
                    cur_depth = jump_depth; // jump to closest bn/parent
                }
            }
        }

        /* // if partial embedding from u in current depth cannot be extended (no complete subgraphs)
        if (partial_embedding_count[cur_depth] == embedding_cnt && cur_depth != 0) {
            int jump_depth = jump[cur_depth];
            for (int i=cur_depth; i>= jump_depth; --i) {
                ui v = embedding[order[i]];
                //reverse_embedding[v] = 0;
                visited_vertices[v] = false; // unvisit and unmap previous query vertex's candidate
            }
            //visited_vertices[ embedding[ order[jump_depth] ] ] = false;
            cur_depth = jump_depth; // jump to closest bn/parent
            continue;
        } */

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else {
            ui v = embedding[order[cur_depth]];
            visited_vertices[v] = false;
        }
    }


    // Release the buffer.
    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    delete[] childs;
    delete[] jump;
    return embedding_cnt;
    //delete[] partial_embedding_count;
    //delete[] reverse_embedding;
}


void
EvaluateQuery::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui *candidates_count, ui *&idx,
                              ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                              ui **&valid_candidate_idx, bool *&visited_vertices) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];

    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID cur_vertex = i;
        ui cur_candidate_num = candidates_count[cur_vertex];

        if (cur_candidate_num > max_candidates_num) {
            max_candidates_num = cur_candidate_num;
        }
    }

    idx = new ui[query_vertices_num];
    idx_count = new ui[query_vertices_num];
    embedding = new ui[query_vertices_num];
    idx_embedding = new ui[query_vertices_num];
    visited_vertices = new bool[data_vertices_num];
    temp_buffer = new ui[max_candidates_num];
    valid_candidate_idx = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        valid_candidate_idx[i] = new ui[max_candidates_num];
    }

    std::fill(visited_vertices, visited_vertices + data_vertices_num, false);
}

void EvaluateQuery::generateValidCandidateIndex(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                                ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                                bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot,
                                                ui **candidates) {
    VertexID u = order[depth];
    VertexID pivot_vertex = pivot[depth];
    ui idx_id = idx_embedding[pivot_vertex];
    Edges &edge = *edge_matrix[pivot_vertex][u];
    ui count = edge.offset_[idx_id + 1] - edge.offset_[idx_id];
    ui *candidate_idx = edge.edge_ + edge.offset_[idx_id];

    ui valid_candidate_index_count = 0;

    if (bn_cnt[depth] == 0) {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];
#ifdef JR
            valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
#else
            if (!visited_vertices[temp_v])
                valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
#endif
        }
    } else {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];
#ifdef JR
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j) {
                VertexID u_bn = bn[depth][j];
                VertexID u_bn_v = embedding[u_bn];

                if (!data_graph->checkEdgeExistence(temp_v, u_bn_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid)
                valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
#else
            if (!visited_vertices[temp_v]) {
                bool valid = true;

                for (ui j = 0; j < bn_cnt[depth]; ++j) {
                    VertexID u_bn = bn[depth][j];
                    VertexID u_bn_v = embedding[u_bn];

                    if (!data_graph->checkEdgeExistence(temp_v, u_bn_v)) {
                        valid = false;
                        break;
                    }
                }

                if (valid)
                    valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
            }
#endif
        }
    }

    idx_count[depth] = valid_candidate_index_count;
}

void EvaluateQuery::generateValidCandidateIndexJR(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                                ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                                ui **bn, ui *bn_cnt, ui *order, ui *pivot, ui **candidates) {
    VertexID u = order[depth];
    VertexID pivot_vertex = pivot[depth];
    ui idx_id = idx_embedding[pivot_vertex];
    Edges &edge = *edge_matrix[pivot_vertex][u];
    ui count = edge.offset_[idx_id + 1] - edge.offset_[idx_id];
    ui *candidate_idx = edge.edge_ + edge.offset_[idx_id];

    ui valid_candidate_index_count = 0;

    if (bn_cnt[depth] == 0) {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            //if (!visited_vertices[temp_v])
            valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
        }
    } else {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            /* if (!visited_vertices[temp_v]) {
            } */
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j) {
                VertexID u_bn = bn[depth][j];
                VertexID u_bn_v = embedding[u_bn];

                if (!data_graph->checkEdgeExistence(temp_v, u_bn_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid)
                valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
        }
    }

    idx_count[depth] = valid_candidate_index_count;
}


void EvaluateQuery::releaseBuffer(ui query_vertices_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn,
                                  ui *bn_count) {
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] idx_embedding;
    delete[] visited_vertices;
    delete[] bn_count;
    delete[] temp_buffer;
    for (ui i = 0; i < query_vertices_num; ++i) {
        delete[] valid_candidate_idx[i];
        delete[] bn[i];
    }

    delete[] valid_candidate_idx;
    delete[] bn;
}

size_t
EvaluateQuery::LFTJ(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count,
                    ui *order, size_t output_limit_num, size_t &call_count) {

#ifdef DISTRIBUTION
    distribution_count_ = new size_t[data_graph->getVerticesCount()];
    memset(distribution_count_, 0, data_graph->getVerticesCount() * sizeof(size_t));
    size_t* begin_count = new size_t[query_graph->getVerticesCount()];
    memset(begin_count, 0, query_graph->getVerticesCount() * sizeof(size_t));
#endif

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET_OTHER
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

#ifdef SPECTRUM
    exit_ = false;
#endif

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            if (visited_vertices[v]) {
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET_OTHER
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

#ifdef DISTRIBUTION
            begin_count[cur_depth] = embedding_cnt;
            // printf("Cur Depth: %d, v: %u, begin: %zu\n", cur_depth, v, embedding_cnt);
#endif

#ifdef ENABLE_FAILING_SET_OTHER
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;

#ifdef DISTRIBUTION
                distribution_count_[v] += 1;
#endif

#ifdef ENABLE_FAILING_SET_OTHER
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                if (embedding_cnt >= output_limit_num) {
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

#ifdef ENABLE_FAILING_SET_OTHER
                if (idx_count[cur_depth] == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                } else {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

#ifdef SPECTRUM
        if (exit_) {
            goto EXIT;
        }
#endif

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET_OTHER
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;

#ifdef DISTRIBUTION
            distribution_count_[embedding[u]] += embedding_cnt - begin_count[cur_depth];
            // printf("Cur Depth: %d, v: %u, begin: %zu, end: %zu\n", cur_depth, embedding[u], begin_count[cur_depth], embedding_cnt);
#endif
        }
    }

    // Release the buffer.

#ifdef DISTRIBUTION
    if (embedding_cnt >= output_limit_num) {
        for (int i = 0; i < max_depth - 1; ++i) {
            ui v = embedding[order[i]];
            distribution_count_[v] += embedding_cnt - begin_count[i];
        }
    }
    delete[] begin_count;
#endif

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i) {
        for (ui j = 0; j < max_depth; ++j) {
//            delete qfliter_bsr_graph_[i][j];
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order,
                                                ui *&temp_buffer) {
    VertexID u = order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;

#if ENABLE_QFLITER == 1
    BSRGraph &bsr_graph = *qfliter_bsr_graph_[previous_bn][u];
    BSRSet &bsr_set = bsr_graph.bsrs[previous_index_id];

    if (bsr_set.size_ != 0){
        offline_bsr_trans_uint(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                               (int *) valid_candidate_index[depth]);
        // represent bsr size
        valid_candidates_count = bsr_set.size_;
    }

    if (bn_cnt[depth] > 0) {
        if (temp_bsr_base1_ == nullptr) { temp_bsr_base1_ = new int[1024 * 1024]; }
        if (temp_bsr_state1_ == nullptr) { temp_bsr_state1_ = new int[1024 * 1024]; }
        if (temp_bsr_base2_ == nullptr) { temp_bsr_base2_ = new int[1024 * 1024]; }
        if (temp_bsr_state2_ == nullptr) { temp_bsr_state2_ = new int[1024 * 1024]; }
        int *res_base_ = temp_bsr_base1_;
        int *res_state_ = temp_bsr_state1_;
        int *input_base_ = temp_bsr_base2_;
        int *input_state_ = temp_bsr_state2_;

        memcpy(input_base_, bsr_set.base_, sizeof(int) * bsr_set.size_);
        memcpy(input_state_, bsr_set.states_, sizeof(int) * bsr_set.size_);

        for (ui i = 1; i < bn_cnt[depth]; ++i) {
            VertexID current_bn = bn[depth][i];
            ui current_index_id = idx_embedding[current_bn];
            BSRGraph &cur_bsr_graph = *qfliter_bsr_graph_[current_bn][u];
            BSRSet &cur_bsr_set = cur_bsr_graph.bsrs[current_index_id];

            if (valid_candidates_count == 0 || cur_bsr_set.size_ == 0) {
                valid_candidates_count = 0;
                break;
            }

            valid_candidates_count = intersect_qfilter_bsr_b4_v2(cur_bsr_set.base_, cur_bsr_set.states_,
                                                                 cur_bsr_set.size_,
                                                                 input_base_, input_state_, valid_candidates_count,
                                                                 res_base_, res_state_);

            swap(res_base_, input_base_);
            swap(res_state_, input_state_);
        }

        if (valid_candidates_count != 0) {
            valid_candidates_count = offline_bsr_trans_uint(input_base_, input_state_, valid_candidates_count,
                                                            (int *) valid_candidate_index[depth]);
        }
    }

    idx_count[depth] = valid_candidates_count;

    // Debugging.
#ifdef YCHE_DEBUG
    Edges &previous_edge = *edge_matrix[previous_bn][u];

    auto gt_valid_candidates_count =
            previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];
    ui *gt_valid_candidate_index = new ui[1024 * 1024];
    memcpy(gt_valid_candidate_index, previous_candidates, gt_valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i) {
        VertexID current_bn = bn[depth][i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
                current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  gt_valid_candidate_index, gt_valid_candidates_count, temp_buffer,
                                                  temp_count);

        std::swap(temp_buffer, gt_valid_candidate_index);
        gt_valid_candidates_count = temp_count;
    }
    assert(valid_candidates_count == gt_valid_candidates_count);

    cout << "Ret, Level:" << bn_cnt[depth] << ", BSR:"
         << pretty_print_array(valid_candidate_index[depth], valid_candidates_count)
         << "; GT: " << pretty_print_array(gt_valid_candidate_index, gt_valid_candidates_count) << "\n";

    for (auto i = 0; i < valid_candidates_count; i++) {
        assert(gt_valid_candidate_index[i] == valid_candidate_index[depth][i]);
    }
    delete[] gt_valid_candidate_index;
#endif
#else
    Edges& previous_edge = *edge_matrix[previous_bn][u];

    valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui* previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i) {
        VertexID current_bn = bn[depth][i];
        Edges& current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui* current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                        temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
#endif
}

size_t EvaluateQuery::exploreGraphQLStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                          ui *candidates_count, ui *order,
                                          size_t output_limit_num, size_t &call_count, TimeOutException* timeout) {
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

#ifdef JR
    // recompute fn/child candidate set if bn/parent's candidate has changed
    int* childs = new int[max_depth]; // fn and child of any query vertex
    int redo = 0; // flag of query vertices that require recomputation of candidate set
    int* jump = new int[max_depth]; // closest bn/parent index to jump to if enumeration failed at any query vertex
    std::vector<int> order_idx(max_depth);
#endif

    // Generate the bn.
    ui **bn;
    ui *bn_count;

    bn = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        bn[i] = new ui[max_depth];
#ifdef JR
        jump[i] = 0;
        childs[i] = 0;
        order_idx[ order[i] ] = i;
#endif
    }

    bn_count = new ui[max_depth];
    std::fill(bn_count, bn_count + max_depth, 0);

    std::vector<bool> visited_query_vertices(max_depth, false);
    visited_query_vertices[start_vertex] = true;
    for (ui i = 1; i < max_depth; ++i) {
        VertexID cur_vertex = order[i];
        ui nbr_cnt;
        const VertexID *nbrs = query_graph->getVertexNeighbors(cur_vertex, nbr_cnt);

#ifdef JR
        int max_idx = 0;
#endif
        for (ui j = 0; j < nbr_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_query_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
#ifdef JR
                int nbr_idx = order_idx[nbr];
                childs[ nbr_idx ] += (1 << i);
                max_idx = std::max(max_idx, nbr_idx);
#endif
            }
        }
#ifdef JR
        jump[i] = max_idx;
#endif

        visited_query_vertices[cur_vertex] = true;
    }

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    for (ui i = 0; i < max_depth; ++i) {
        VertexID cur_vertex = order[i];
        ui max_candidate_count = candidates_count[cur_vertex];
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;
#ifdef JR
            if (visited_vertices[v]) continue;
            // REDO
            redo |= childs[cur_depth];
#endif

            embedding[u] = v;
            visited_vertices[v] = true;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                cur_depth += 1;
                idx[cur_depth] = 0;
#ifdef JR
                if (redo & (1 << cur_depth) ) {
                    redo -= (1 << cur_depth);
                    call_count += 1;
                    generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                            visited_vertices, bn, bn_count, order, candidates, candidates_count);
                }

                // JUMP
                if (idx_count[cur_depth] == 0) {
                    int jump_depth = jump[cur_depth];
                    for (int i=cur_depth-1; i>= jump_depth; --i) {
                        ui v = embedding[order[i]];
                        visited_vertices[v] = false; // unvisit and unmap previous query vertex's candidate
                    }
                    cur_depth = jump_depth; // jump to closest bn/parent
                }
#else
                call_count += 1;
                generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, candidates, candidates_count);
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

    // Release the buffer.
    EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

#ifdef JR
    delete[] childs;
    delete[] jump;
#endif

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidates(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                            ui **valid_candidate, bool *visited_vertices, ui **bn, ui *bn_cnt,
                                            ui *order, ui **candidates, ui *candidates_count) {
    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i) {
        VertexID v = candidates[u][i];

        
#ifdef JR
        bool valid = true;

        for (ui j = 0; j < bn_cnt[depth]; ++j) {
            VertexID u_nbr = bn[depth][j];
            VertexID u_nbr_v = embedding[u_nbr];

            if (!data_graph->checkEdgeExistence(v, u_nbr_v)) {
                valid = false;
                break;
            }
        }

        if (valid) {
            valid_candidate[depth][idx_count[depth]++] = v;
        }
#else
        if (!visited_vertices[v]) {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j) {
                VertexID u_nbr = bn[depth][j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
#endif
    }
}

size_t EvaluateQuery::exploreQuickSIStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                          ui *candidates_count, ui *order,
                                          ui *pivot, size_t output_limit_num, size_t &call_count) {
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    // Generate the bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, pivot, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num) {
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(query_graph, data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, pivot);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

    // Release the buffer.
    EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidates(const Graph *query_graph, const Graph *data_graph, ui depth, ui *embedding,
                                            ui *idx_count, ui **valid_candidate, bool *visited_vertices, ui **bn,
                                            ui *bn_cnt,
                                            ui *order, ui *pivot) {
    VertexID u = order[depth];
    LabelID u_label = query_graph->getVertexLabel(u);
    ui u_degree = query_graph->getVertexDegree(u);

    idx_count[depth] = 0;

    VertexID p = embedding[pivot[depth]];
    ui nbr_cnt;
    const VertexID *nbrs = data_graph->getVertexNeighbors(p, nbr_cnt);

    for (ui i = 0; i < nbr_cnt; ++i) {
        VertexID v = nbrs[i];

        if (!visited_vertices[v] && u_label == data_graph->getVertexLabel(v) &&
            u_degree <= data_graph->getVertexDegree(v)) {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j) {
                VertexID u_nbr = bn[depth][j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}

size_t EvaluateQuery::exploreDPisoStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                        Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                        ui **weight_array, ui *order, size_t output_limit_num,
                                        size_t &call_count, TimeOutException* timeout) {
    int max_depth = query_graph->getVerticesCount();

    ui *extendable = new ui[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        extendable[i] = tree[i].bn_count_;
    }

    // Generate backward neighbors.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    VertexID start_vertex = order[0];
    std::vector<dpiso_min_pq> vec_rank_queue;

    for (ui i = 0; i < candidates_count[start_vertex]; ++i) {
        VertexID v = candidates[start_vertex][i];

        embedding[start_vertex] = v;
        idx_embedding[start_vertex] = i;
        visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
        reverse_embedding[v] = start_vertex;
#endif
        vec_rank_queue.emplace_back(dpiso_min_pq(extendable_vertex_compare));
        updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer, weight_array,
                               tree, start_vertex, extendable,
                               vec_rank_queue, query_graph);

        VertexID u = vec_rank_queue.back().top().first.first;
        vec_rank_queue.back().pop();

#ifdef ENABLE_FAILING_SET
        if (idx_count[u] == 0) {
            vec_failing_set[cur_depth] = ancestors[u];
        } else {
            vec_failing_set[cur_depth].reset();
        }
#endif

        call_count += 1;
        cur_depth += 1;
        order[cur_depth] = u;
        idx[u] = 0;
        while (cur_depth > 0) {
            while (idx[u] < idx_count[u]) {
                ui valid_idx = valid_candidate_idx[u][idx[u]];
                v = candidates[u][valid_idx];

                if (visited_vertices[v]) {
                    idx[u] += 1;
#ifdef ENABLE_FAILING_SET
                    vec_failing_set[cur_depth] = ancestors[u];
                    vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                    continue;
                }

                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_vertices[v] = true;
                idx[u] += 1;

#ifdef ENABLE_FAILING_SET
                reverse_embedding[v] = u;
#endif

                if (cur_depth == max_depth - 1) {
                    embedding_cnt += 1;
                    visited_vertices[v] = false;
#ifdef ENABLE_FAILING_SET
                    reverse_embedding.erase(embedding[u]);
                    vec_failing_set[cur_depth].set();
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

#endif
                    if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                        goto EXIT;
                    }
                } else {
                    cur_depth += 1;

                    call_count += 1;
                    vec_rank_queue.emplace_back(vec_rank_queue.back());
                    updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer,
                                           weight_array, tree, u, extendable,
                                           vec_rank_queue, query_graph);

                    u = vec_rank_queue.back().top().first.first;
                    vec_rank_queue.back().pop();
                    idx[u] = 0;
                    order[cur_depth] = u;

#ifdef ENABLE_FAILING_SET
                    if (idx_count[u] == 0) {
                        vec_failing_set[cur_depth - 1] = ancestors[u];
                    } else {
                        vec_failing_set[cur_depth - 1].reset();
                    }
#endif
                }
            }

            cur_depth -= 1;
            vec_rank_queue.pop_back();
            u = order[cur_depth];
            visited_vertices[embedding[u]] = false;
            restoreExtendableVertex(tree, u, extendable);
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[u] = idx_count[u];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            if (timeout->is_throwable)
                goto EXIT;
        }
    }

    // Release the buffer.
    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

void EvaluateQuery::updateExtendableVertex(ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                           Edges ***edge_matrix, ui *&temp_buffer, ui **weight_array,
                                           TreeNode *tree, VertexID mapped_vertex, ui *extendable,
                                           std::vector<dpiso_min_pq> &vec_rank_queue, const Graph *query_graph) {
    TreeNode &node = tree[mapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i) {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0) {
            generateValidCandidateIndex(u, idx_embedding, idx_count, valid_candidate_index[u], edge_matrix, tree[u].bn_,
                                        tree[u].bn_count_, temp_buffer);

            ui weight = 0;
            for (ui j = 0; j < idx_count[u]; ++j) {
                ui idx = valid_candidate_index[u][j];
                weight += weight_array[u][idx];
            }
            vec_rank_queue.back().emplace(std::make_pair(std::make_pair(u, query_graph->getVertexDegree(u)), weight));
        }
    }
}

void EvaluateQuery::restoreExtendableVertex(TreeNode *tree, VertexID unmapped_vertex, ui *extendable) {
    TreeNode &node = tree[unmapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i) {
        VertexID u = node.fn_[i];
        extendable[u] += 1;
    }
}

void
EvaluateQuery::generateValidCandidateIndex(VertexID u, ui *idx_embedding, ui *idx_count, ui *&valid_candidate_index,
                                           Edges ***edge_matrix, ui *bn, ui bn_cnt, ui *&temp_buffer) {
    VertexID previous_bn = bn[0];
    Edges &previous_edge = *edge_matrix[previous_bn][u];
    ui previous_index_id = idx_embedding[previous_bn];

    ui previous_candidates_count =
            previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    ui valid_candidates_count = 0;
    for (ui i = 0; i < previous_candidates_count; ++i) {
        valid_candidate_index[valid_candidates_count++] = previous_candidates[i];
    }

    ui temp_count;
    for (ui i = 1; i < bn_cnt; ++i) {
        VertexID current_bn = bn[i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
                current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index,
                                                  valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index);
        valid_candidates_count = temp_count;
    }

    idx_count[u] = valid_candidates_count;
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, TreeNode *tree, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        TreeNode &u_node = tree[u];
        ancestors[u].set(u);
        for (ui j = 0; j < u_node.bn_count_; ++j) {
            VertexID u_bn = u_node.bn_[j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}

size_t EvaluateQuery::exploreDPisoRecursiveStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                                 Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                                 ui **weight_array, ui *order, size_t output_limit_num,
                                                 size_t &call_count) {
    int max_depth = query_graph->getVerticesCount();

    ui *extendable = new ui[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        extendable[i] = tree[i].bn_count_;
    }

    // Generate backward neighbors.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    // Evaluate the query.
    size_t embedding_cnt = 0;

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);

    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    VertexID start_vertex = order[0];

    for (ui i = 0; i < candidates_count[start_vertex]; ++i) {
        VertexID v = candidates[start_vertex][i];
        embedding[start_vertex] = v;
        idx_embedding[start_vertex] = i;
        visited_vertices[v] = true;
        reverse_embedding[v] = start_vertex;
        call_count += 1;

        exploreDPisoBacktrack(max_depth, 1, start_vertex, tree, idx_embedding, embedding, reverse_embedding,
                              visited_vertices, idx_count, valid_candidate_idx, edge_matrix,
                              ancestors, dpiso_min_pq(extendable_vertex_compare), weight_array, temp_buffer, extendable,
                              candidates, embedding_cnt,
                              call_count, nullptr);

        visited_vertices[v] = false;
        reverse_embedding.erase(v);
    }

    // Release the buffer.
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>
EvaluateQuery::exploreDPisoBacktrack(ui max_depth, ui depth, VertexID mapped_vertex, TreeNode *tree, ui *idx_embedding,
                                     ui *embedding, std::unordered_map<VertexID, VertexID> &reverse_embedding,
                                     bool *visited_vertices, ui *idx_count, ui **valid_candidate_index,
                                     Edges ***edge_matrix,
                                     std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors,
                                     dpiso_min_pq rank_queue, ui **weight_array, ui *&temp_buffer, ui *extendable,
                                     ui **candidates, size_t &embedding_count, size_t &call_count,
                                     const Graph *query_graph) {
    // Compute extendable vertex.
    TreeNode &node = tree[mapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i) {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0) {
            generateValidCandidateIndex(u, idx_embedding, idx_count, valid_candidate_index[u], edge_matrix, tree[u].bn_,
                                        tree[u].bn_count_, temp_buffer);

            ui weight = 0;
            for (ui j = 0; j < idx_count[u]; ++j) {
                ui idx = valid_candidate_index[u][j];
                weight += weight_array[u][idx];
            }
            rank_queue.emplace(std::make_pair(std::make_pair(u, query_graph->getVertexDegree(u)), weight));
        }
    }

    VertexID u = rank_queue.top().first.first;
    rank_queue.pop();

    if (idx_count[u] == 0) {
        restoreExtendableVertex(tree, mapped_vertex, extendable);
        return ancestors[u];
    } else {
        std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> current_fs;
        std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> child_fs;

        for (ui i = 0; i < idx_count[u]; ++i) {
            ui valid_index = valid_candidate_index[u][i];
            VertexID v = candidates[u][valid_index];

            if (!visited_vertices[v]) {
                embedding[u] = v;
                idx_embedding[u] = valid_index;
                visited_vertices[v] = true;
                reverse_embedding[v] = u;
                if (depth != max_depth - 1) {
                    call_count += 1;
                    child_fs = exploreDPisoBacktrack(max_depth, depth + 1, u, tree, idx_embedding, embedding,
                                                     reverse_embedding, visited_vertices, idx_count,
                                                     valid_candidate_index, edge_matrix,
                                                     ancestors, rank_queue, weight_array, temp_buffer, extendable,
                                                     candidates, embedding_count,
                                                     call_count, query_graph);
                } else {
                    embedding_count += 1;
                    child_fs.set();
                }
                visited_vertices[v] = false;
                reverse_embedding.erase(v);

                if (!child_fs.test(u)) {
                    current_fs = child_fs;
                    break;
                }
            } else {
                child_fs.reset();
                child_fs |= ancestors[u];
                child_fs |= ancestors[reverse_embedding[v]];
            }

            current_fs |= child_fs;
        }

        restoreExtendableVertex(tree, mapped_vertex, extendable);
        return current_fs;
    }
}

size_t
EvaluateQuery::exploreCECIStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout) {

    int max_depth = query_graph->getVerticesCount();
#ifdef JR
    // recompute fn/child candidate set if bn/parent's candidate has changed
    int* childs = new int[max_depth]; // fn and child of any query vertex
    int redo = 0; // flag of query vertices that require recomputation of candidate set
    int* jump = new int[max_depth]; // closest bn/parent index to jump to if enumeration failed at any query vertex
    std::vector<int> order_idx(max_depth);
#endif

    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (int i = 0; i < max_depth; ++i) {
        if (candidates_count[i] > max_valid_candidates_count) {
            max_valid_candidates_count = candidates_count[i];
        }
#ifdef JR
        jump[i] = 0;
        childs[i] = 0;
        order_idx[ order[i] ] = i;
#endif
    }
    // Allocate the memory buffer.
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (int i = 0; i < max_depth; ++i) {
        valid_candidates[i] = new ui[max_valid_candidates_count];

#ifdef JR
        // JR initialization
        VertexID u = order[i];
        int max_idx = 0;
        ui deg = 0;

        const ui* nbrs = query_graph->getVertexNeighbors(u, deg);
        for (int j=0; j<deg; ++j) {
            int u_b_idx = order_idx[ nbrs[j] ];
            if (u_b_idx > i)
                continue;
            childs[ u_b_idx ] += (1 << i);
            max_idx = std::max(max_idx, u_b_idx);
        }

        jump[i] = max_idx;
#endif

    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

#ifdef JR
            // REDO
            redo |= childs[cur_depth];
#endif
            embedding[u] = v;
            visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif
            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                cur_depth += 1;
                idx[cur_depth] = 0;

#ifdef JR
                if (redo & (1 << cur_depth) ) {
                    redo -= (1 << cur_depth);

                    call_count += 1;
                    idx_count[cur_depth] = 0;
                    generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                            TE_Candidates,
                                            NTE_Candidates);
                }

                // JUMP
                if (idx_count[cur_depth] == 0) {
                    int jump_depth = jump[cur_depth];
                    for (int i=cur_depth-1; i>= jump_depth; --i) {
                        ui v = embedding[order[i]];
                        visited_vertices[v] = false; // unvisit and unmap previous query vertex's candidate
                    }
                    cur_depth = jump_depth; // jump to closest bn/parent
                }
#else
                call_count += 1;
                generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                        TE_Candidates,
                                        NTE_Candidates);
#endif

#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                } else {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.
    EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (int i = 0; i < max_depth; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;

#ifdef JR
    delete[] childs;
    delete[] jump;
#endif
    return embedding_cnt;
}

// l2Match with JR
size_t
EvaluateQuery::exploreL2Match(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout) {

    int max_depth = query_graph->getVerticesCount();

    // JR variables
    // recompute fn/child candidate set if bn/parent's candidate has changed
    int* childs = new int[max_depth]; // fn and child of any query vertex
    int redo = 0; // flag of query vertices that require recomputation of candidate set
    int* jump = new int[max_depth]; // closest bn/parent index to jump to if enumeration failed at any query vertex
    std::vector<int> order_idx(max_depth);

    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (int i = 0; i < max_depth; ++i) {
        if (candidates_count[i] > max_valid_candidates_count) {
            max_valid_candidates_count = candidates_count[i];
        }
        jump[i] = 0;
        childs[i] = 0;
        order_idx[ order[i] ] = i;
    }
    // Allocate the memory buffer.
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (int i = 0; i < max_depth; ++i) {
        valid_candidates[i] = new ui[max_valid_candidates_count];
        
        // JR initialization
        VertexID u = order[i];
        int max_idx = 0;
        ui deg = 0;

        const ui* nbrs = query_graph->getVertexNeighbors(u, deg);
        for (int j=0; j<deg; ++j) {
            int u_b_idx = order_idx[ nbrs[j] ];
            if (u_b_idx > i)
                continue;
            childs[ u_b_idx ] += (1 << i);
            max_idx = std::max(max_idx, u_b_idx);
        }

        jump[i] = max_idx;
    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {            
                continue;
            }

            // REDO
            redo |= childs[cur_depth];

            embedding[u] = v;
            visited_vertices[v] = true;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;

                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                cur_depth += 1;
                idx[cur_depth] = 0;
                
                if (redo & (1 << cur_depth) ) {
                    redo -= (1 << cur_depth);

                    call_count += 1;
                    idx_count[cur_depth] = 0;
                    generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                            TE_Candidates,
                                            NTE_Candidates);
                }

                // JUMP
                if (idx_count[cur_depth] == 0) {
                    int jump_depth = jump[cur_depth];
                    for (int i=cur_depth-1; i>= jump_depth; --i) {
                        ui v = embedding[order[i]];
                        visited_vertices[v] = false; // unvisit and unmap previous query vertex's candidate
                    }
                    cur_depth = jump_depth; // jump to closest bn/parent
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else {
            VertexID u = order[cur_depth];
          
            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.
    EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (int i = 0; i < max_depth; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;
    delete[] childs;
    delete[] jump;

    return embedding_cnt;
}

// l2Match w/o TE_Candidates & with Proximity Ordering
size_t
EvaluateQuery::exploreL2MatchB(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout) {

    int max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (int i = 0; i < max_depth; ++i) {
        if (candidates_count[i] > max_valid_candidates_count) {
            max_valid_candidates_count = candidates_count[i];
        }
    }
    // Allocate the memory buffer.
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (int i = 0; i < max_depth; ++i) {
        valid_candidates[i] = new ui[max_valid_candidates_count];
    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                cur_depth += 1;
                call_count += 1;
                idx[cur_depth] = 0;
                generateValidCandidatesL2MatchB(cur_depth, embedding, idx_count, valid_candidates, order, tree, NTE_Candidates);

                
#ifdef ENABLE_FAILING_SET                
                if (idx_count[cur_depth] == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                } else {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif   
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else {
            VertexID u = order[cur_depth];

#ifdef ENABLE_FAILING_SET            
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif

            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.
    EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (int i = 0; i < max_depth; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;

    return embedding_cnt;
}

// l2Match w/o JR
size_t
EvaluateQuery::exploreL2MatchC(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, size_t &output_limit_num, size_t &call_count, TimeOutException* timeout) {

    int max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (int i = 0; i < max_depth; ++i) {
        if (candidates_count[i] > max_valid_candidates_count) {
            max_valid_candidates_count = candidates_count[i];
        }
    }
    // Allocate the memory buffer.
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (int i = 0; i < max_depth; ++i) {
        valid_candidates[i] = new ui[max_valid_candidates_count];
    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {
#ifdef ENABLE_FAILING_SET                
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif                
                continue;
            }

            embedding[u] = v;
            visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;

#ifdef ENABLE_FAILING_SET                
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                        TE_Candidates,
                                        NTE_Candidates);

#ifdef ENABLE_FAILING_SET                                        
                if (idx_count[cur_depth] == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                } else {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else {
            VertexID u = order[cur_depth];

#ifdef ENABLE_FAILING_SET            
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif            
            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.
    EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (int i = 0; i < max_depth; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;

    return embedding_cnt;
}




void EvaluateQuery::generateValidCandidates(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            ui *&temp_buffer, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {

    VertexID u = order[depth];
    TreeNode &u_node = tree[u];
    idx_count[depth] = 0;
    ui valid_candidates_count = 0;
    {
        VertexID u_p = tree[u].parent_;
        VertexID v_p = embedding[u_p];

        auto iter = TE_Candidates[u].find(v_p);
        if (iter == TE_Candidates[u].end() || iter->second.empty()) {
            return;
        }

        valid_candidates_count = iter->second.size();
        VertexID *v_p_nbrs = iter->second.data();

        for (ui i = 0; i < valid_candidates_count; ++i) {
            valid_candidates[depth][i] = v_p_nbrs[i];
        }
    }
    ui temp_count;
    for (ui i = 0; i < tree[u].bn_count_; ++i) {
        VertexID u_p = tree[u].bn_[i];
        VertexID v_p = embedding[u_p];

        auto iter = NTE_Candidates[u][u_p].find(v_p);
        if (iter == NTE_Candidates[u][u_p].end() || iter->second.empty()) {
            return;
        }

        ui current_candidates_count = iter->second.size();
        ui *current_candidates = iter->second.data();

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  valid_candidates[depth], valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidates[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}

// l2Match w/o TE_Candidates & with Proximity Ordering
void EvaluateQuery::generateValidCandidatesL2MatchB(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order, TreeNode *tree,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {

    VertexID u = order[depth];
    TreeNode &u_node = tree[u];
    idx_count[depth] = 0;
    ui count = 0;
    {
        VertexID u_p = u_node.parent_;
        VertexID v_p = embedding[u_p];

        auto iter = NTE_Candidates[u][u_p].find(v_p);
        if (iter == NTE_Candidates[u][u_p].end() || iter->second.empty()) {
            return;
        }

        ui v_p_nbrs_count = iter->second.size();
        VertexID *v_p_nbrs = iter->second.data();

        for (ui i = 0; i < v_p_nbrs_count; ++i) {
            ui temp_v = v_p_nbrs[i];
            valid_candidates[depth][count++] = temp_v;
        }

        if (u_node.bn_count_ == 0) {
            idx_count[depth] = count;
            return;
        }
    }
    
    ui valid_candidates_count = 0;
    for (ui i = 0; i < count; ++i) {
        VertexID temp_v = valid_candidates[depth][i];

        bool valid = true;

        for (ui j = 0; j < u_node.bn_count_; ++j) {
            VertexID u_bn = u_node.bn_[j];
            VertexID u_bn_v = embedding[u_bn];

            auto iter = NTE_Candidates[u][u_bn].find(u_bn_v);
            if (iter == NTE_Candidates[u][u_bn].end() || iter->second.empty()) {
                valid = false;
                break;
            }
            if (!binary_search(iter->second.begin(), iter->second.end(), temp_v)) {
                valid = false;
                break;
            }
        }

        if (valid)
            valid_candidates[depth][valid_candidates_count++] = temp_v;


    }

    idx_count[depth] = valid_candidates_count;
}


void EvaluateQuery::generateValidCandidatesL2MatchC(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates,
                                            ui *order, TreeNode *tree, bool *visited_vertices,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {

    VertexID u = order[depth];
    TreeNode &u_node = tree[u];
    idx_count[depth] = 0;
    ui count = 0;
    {
        VertexID u_p = u_node.parent_;
        VertexID v_p = embedding[u_p];

        auto iter = NTE_Candidates[u][u_p].find(v_p);
        if (iter == NTE_Candidates[u][u_p].end() || iter->second.empty()) {
            return;
        }

        ui v_p_nbrs_count = iter->second.size();
        VertexID *v_p_nbrs = iter->second.data();

        for (ui i = 0; i < v_p_nbrs_count; ++i) {
            ui temp_v = v_p_nbrs[i];
            if (!visited_vertices[temp_v])
                valid_candidates[depth][count++] = temp_v;
        }

        if (u_node.bn_count_ == 0) {
            idx_count[depth] = count;
            return;
        }
    }
    
    ui valid_candidates_count = 0;
    for (ui i = 0; i < count; ++i) {
        VertexID temp_v = valid_candidates[depth][i];

        if (!visited_vertices[temp_v]) {
            bool valid = true;

            for (ui j = 0; j < u_node.bn_count_; ++j) {
                VertexID u_bn = u_node.bn_[j];
                VertexID u_bn_v = embedding[u_bn];

                auto iter = NTE_Candidates[u][u_bn].find(u_bn_v);
                if (iter == NTE_Candidates[u][u_bn].end() || iter->second.empty()) {
                    valid = false;
                    break;
                }
                if (!binary_search(iter->second.begin(), iter->second.end(), temp_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid)
                valid_candidates[depth][valid_candidates_count++] = temp_v;

        }


    }

    idx_count[depth] = valid_candidates_count;
}


void EvaluateQuery::computeAncestor(const Graph *query_graph, ui **bn, ui *bn_cnt, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < bn_cnt[i]; ++j) {
            VertexID u_bn = bn[i][j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < i; ++j) {
            VertexID u_bn = order[j];
            if (query_graph->checkEdgeExistence(u, u_bn)) {
                ancestors[u] |= ancestors[u_bn];
            }
        }
    }
}


// l2Match with JR without CECI index
size_t
EvaluateQuery::exploreL2MatchNew(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, TreeNode *tree, ui **candidates,
                            ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count, TimeOutException* timeout) {
    // Generate the bn.
    int max_depth = query_graph->getVerticesCount();
    // recompute fn/child candidate set if bn/parent's candidate has changed
    int* childs = new int[max_depth]; // fn and child of any query vertex
    int redo = 0; // flag of query vertices that require recomputation of candidate set
    int* jump = new int[max_depth]; // closest bn/parent index to jump to if enumeration failed at any query vertex
    std::vector<int> order_idx(max_depth);

    ui **bn = new ui*[max_depth];
    ui *bn_count = new ui[max_depth];
    pivot = new ui[max_depth];

    for (int i = 0; i < max_depth; ++i) {
        jump[i] = 0;
        childs[i] = 0;
        order_idx[ order[i] ] = i;
    }

    for (int i = 0; i < max_depth; ++i) {
        // JR initialization
        VertexID u = order[i];
        pivot[i] = tree[u].parent_;
        int max_idx = 0;
        ui deg = 0;
        bn[i] = new ui[max_depth];

        const ui* nbrs = query_graph->getVertexNeighbors(u, deg);
        for (int j=0; j<deg; ++j) {
            int u_b_idx = order_idx[ nbrs[j] ];
            if (u_b_idx > i)
                continue;
            childs[ u_b_idx ] += (1 << i);
            max_idx = std::max(max_idx, u_b_idx);
            bn[i][ bn_count[i]++ ] = nbrs[j];
        }

        jump[i] = max_idx;
    }

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {
                // find next non-visited valid candidate
                continue;
            }
            // REDO
            redo |= childs[cur_depth];

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num || timeout->is_throwable) {
                    goto EXIT;
                }
            } else {
                cur_depth += 1;
                idx[cur_depth] = 0;
                
                if (redo & (1 << cur_depth)) {
                    redo -= (1 << cur_depth);

                    call_count += 1;
                    idx_count[cur_depth] = 0;
                    generateValidCandidateIndexJRNew(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx, edge_matrix, bn, bn_count,
                                            order, pivot[cur_depth], candidates);
                }

                if (idx_count[cur_depth] == 0) {
                    int jump_depth = jump[cur_depth];
                    for (int i=cur_depth-1; i>= jump_depth; --i) {
                        ui v = embedding[order[i]];
                        visited_vertices[v] = false; // unvisit and unmap previous query vertex's candidate
                    }
                    cur_depth = jump_depth; // jump to closest bn/parent
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0 || timeout->is_throwable)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }


    // Release the buffer.
    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    delete[] childs;
    delete[] jump;

    return embedding_cnt;
}


void EvaluateQuery::generateValidCandidateIndexJRNew(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                                ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                                ui **bn, ui *bn_cnt, ui *order, ui parent, ui **candidates) {
    VertexID u = order[depth];
    ui idx_id = idx_embedding[parent];
    Edges &edge = *edge_matrix[parent][u];
    ui count = edge.offset_[idx_id + 1] - edge.offset_[idx_id];
    ui *candidate_idx = edge.edge_ + edge.offset_[idx_id];

    ui valid_candidate_index_count = 0;

    if (bn_cnt[depth] == 0) {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            //if (!visited_vertices[temp_v])
            valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
        }
    } else {
        for (ui i = 0; i < count; ++i) {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            /* if (!visited_vertices[temp_v]) {
            } */
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j) {
                VertexID u_bn = bn[depth][j];
                VertexID u_bn_v = embedding[u_bn];

                if (!data_graph->checkEdgeExistence(temp_v, u_bn_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid)
                valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
        }
    }

    idx_count[depth] = valid_candidate_index_count;
}


size_t
EvaluateQuery::exploreUllman(const Graph* data_graph, const Graph* query_graph, size_t output_limit_num, TimeOutException* timeout)
{
    ui max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    bool* visited_candidates = new bool[data_vertices_count];
    std::fill(visited_candidates, visited_candidates + data_vertices_count, false);
    ui call_count = 0;
    size_t embedding_count = 0;

    ui** M = new ui*[max_depth];
    for (int i=0; i<max_depth; ++i)
    {
        M[i] = new ui[data_vertices_count];
        for (int j=0; j<data_vertices_count; ++j)
            std::fill(M[i], M[i] + data_vertices_count, 0);
    }

    recurseUllman(visited_candidates, 0, data_graph, query_graph, M, call_count, embedding_count, output_limit_num, timeout);

    delete[] visited_candidates;
    for (int i=0; i<max_depth; ++i)
        delete[] M[i];
    delete[] M;

    return embedding_count;
}

bool
EvaluateQuery::recurseUllman(bool* visited_candidates, ui depth, const Graph* data_graph, const Graph* query_graph, ui** mappings, ui& call_count, size_t& embedding_count, size_t output_limit_num, TimeOutException* timeout)
{
    ui max_depth = query_graph->getVerticesCount();
    if (depth == max_depth)
    {
        embedding_count += 1;
        return true;
    }

    call_count += 1;
    ui data_vertices_count = data_graph->getVerticesCount();
    ui** M_copy = new ui*[max_depth];
    for (int i=0; i<max_depth; ++i)
        M_copy[i] = new ui[data_vertices_count];
    for (int i=0; i<max_depth; ++i)
        std::memcpy(M_copy[i], mappings[i], sizeof(ui) * data_vertices_count);

    pruneUllman(data_graph, query_graph, M_copy);

    bool is_isomorphic = false;
    for (ui col=0; col<data_vertices_count; ++col)
    {
        if (visited_candidates[col])
            continue;

        // set column c in M' to 1 and other columns to 0
        M_copy[depth][col] = 1;
        visited_candidates[col] = true;
        is_isomorphic = recurseUllman(visited_candidates, depth+1, data_graph, query_graph, M_copy, call_count, embedding_count, output_limit_num, timeout);
        visited_candidates[col] = false;
        
        if (embedding_count >= output_limit_num || timeout->is_throwable)
            goto EXIT;
    }

    EXIT:
    return is_isomorphic;
}

void
EvaluateQuery::pruneUllman(const Graph* data_graph, const Graph* query_graph, ui** mappings)
{
    ui max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();

    bool is_changed = false;
    do
    {
        for (int i=0; i<max_depth; ++i)
        {
            ui u_deg = 0;
            const ui* u_nbrs = query_graph->getVertexNeighbors(i, u_deg);

            for (int j=0; j<data_vertices_count; ++j)
            {
                if (mappings[i][j] == 1)
                {
                    ui v_deg = 0;
                    const ui* v_nbrs = data_graph->getVertexNeighbors(j, v_deg);
                    
                    for (ui k=0; k<u_deg; ++k)
                    {
                        ui u_nbr = u_nbrs[k];
                        bool is_valid = false;

                        for (ui l=0; l<v_deg; ++l)
                        {
                            if (mappings[u_nbr][ v_nbrs[l] ] == 1)
                            {
                                is_valid = true;
                                break;
                            }
                        }

                        if (!is_valid)
                        {
                            is_changed = true;
                            mappings[i][j] = 0;
                        }
                    }
                }
            }
        }
    }
    while (is_changed);
}



#endif //SUBGRAPHMATCHING_EVALUATEQUERY_H
