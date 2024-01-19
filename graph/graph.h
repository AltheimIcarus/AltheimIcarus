//
// Created by ssunah on 6/22/18.
//

#ifndef SUBGRAPHMATCHING_GRAPH_H
#define SUBGRAPHMATCHING_GRAPH_H

#include <unordered_map>
#include <vector>
#include <iostream>

#include "configuration/types.h"
#include "configuration/config.h"

/**
 * A graph is stored as the CSR format.
 */

class Graph {
private:
    bool enable_label_offset_;

    ui vertices_count_;
    ui edges_count_;
    ui labels_count_;
    ui max_degree_;
    ui max_label_frequency_;
    double avg_deg_;

    ui* offsets_;
    VertexID * neighbors_;
    VertexID* neighbors_by_labels_;
    LabelID* labels_;
    ui* reverse_index_offsets_;
    ui* reverse_index_;

    int* core_table_;
    ui core_length_;

    std::unordered_map<LabelID, ui> labels_frequency_;

    ui* spatial_index_;
    ui max_spatial_index_;

#if OPTIMIZED_LABELED_GRAPH == 1
    ui* labels_offsets_;
    std::unordered_map<LabelID, ui>* nlf_;
    std::unordered_map<LabelID, std::pair<ui,ui>>* nlo; // nlo[u][label] = (start index, count) for list of adjs with same label id
    // Label Pair Index [L(u)][L(v)] = [u'], v' is accessible via CSR thus not required in the label pair index
public:
    std::vector<std::unordered_map<LabelID, std::vector<VertexID >>> label_pairs;

#endif

private:
    void BuildReverseIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    void BuildNLF();
    void BuildLabelOffset();
#endif

public:
    Graph(const bool enable_label_offset) {
        enable_label_offset_ = enable_label_offset;

        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;
        avg_deg_ = 0.0;
        core_length_ = 0;

        offsets_ = NULL;
        neighbors_ = NULL;
        neighbors_by_labels_ = NULL;
        labels_ = NULL;
        reverse_index_offsets_ = NULL;
        reverse_index_ = NULL;
        core_table_ = NULL;
        labels_frequency_.clear();

        spatial_index_ = NULL;
        max_spatial_index_ = 0;

#if OPTIMIZED_LABELED_GRAPH == 1
        labels_offsets_ = NULL;
        nlf_ = NULL;
        nlo = NULL;
#endif
    }

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] neighbors_by_labels_;
        delete[] labels_;
        delete[] reverse_index_offsets_;
        delete[] reverse_index_;
        delete[] core_table_;
#if OPTIMIZED_LABELED_GRAPH == 1
        delete[] labels_offsets_;
        delete[] spatial_index_;
        delete[] nlf_;
        delete[] nlo;
#endif
    }

public:
    void loadGraphFromFile(const std::string& file_path);
    void loadGraphFromVertexFile(const std::string &vertex_path, const std::string &edge_path);
    void loadGraphFromEdgeFile(const std::string& file_path);
    void loadGraphFromFileBeta(const std::string& file_path);
    void loadGraphFromFileCompressed(const std::string& degree_path, const std::string& edge_path,
                                     const std::string& label_path);
    void storeComparessedGraph(const std::string& degree_path, const std::string& edge_path,
                               const std::string& label_path);
    void printGraphMetaData();
public:
    const double getAverageDegree() const {
        return avg_deg_;
    }
    const ui getLabelsCount() const {
        return labels_count_;
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return edges_count_;
    }

    const ui getGraphMaxDegree() const {
        return max_degree_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const ui getCoreValue(const VertexID id) const {
        return core_table_[id];
    }

    const ui get2CoreSize() const {
        return core_length_;
    }
    const LabelID getVertexLabel(const VertexID id) const {
        return labels_[id];
    }

    const ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }


    const ui * getVerticesByLabel(const LabelID id, ui& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

    const ui getDensity() const {
        return 2 * edges_count_ / (vertices_count_ * (vertices_count_ - 1));
    }

#if OPTIMIZED_LABELED_GRAPH == 1
    // l2Match
    const ui * getNeighborsByLabel(const VertexID id, const LabelID neighbor_label, ui& count) const {
        auto it = nlo[id].find(neighbor_label);
        count = it->second.second;
        return neighbors_by_labels_ + it->second.first;
    }

    const std::unordered_map<LabelID, ui>* getVertexNLF(const VertexID id) const {
        return nlf_ + id;
    }
    // l2Match
    const ui getVertexLabelFrequency(const ui id, const LabelID neighbor_label) const {
        auto it = nlo[id].find(neighbor_label);
        if (it != nlo[id].end())
            return it->second.second;
        return 0;
    }
    // l2Match
    const ui getLabelPairFrequency(const LabelID u_label, const LabelID neighbor_label) const {
        auto it = label_pairs[u_label].find(neighbor_label);
        if (it != label_pairs[u_label].end())
            return it->second.size();
        return 0;
    }
#if OPTIMIZED_LABELED_GRAPH == 1
    const ui getLabelPairFrequency2(const LabelID u_label, const LabelID neighbor_label) const {
        ui count;
        const ui* vertices_by_label = getVerticesByLabel(u_label, count);

        ui freq = 0;
        for (ui i=0; i<count; ++i) {
            ui u = vertices_by_label[i];

            const std::unordered_map<LabelID, ui>* this_vertex_nlf = getVertexNLF(u);
            auto it = this_vertex_nlf->find(neighbor_label);
            if (it != this_vertex_nlf->end())
                freq += 1;
        }
        return freq;
    }
#endif
    // l2Match
    const std::unordered_map<LabelID, std::vector<ui>>::const_iterator getLabelPairIterator(const LabelID u_label, const LabelID neighbor_label) const {
        return label_pairs[u_label].find(neighbor_label);
    }

#endif

    bool checkEdgeExistence(VertexID u, VertexID v, bool search_by_label = false) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;

        const VertexID* neighbors = NULL;
        if (search_by_label) {
            neighbors =  getNeighborsByLabel(v, getVertexLabel(u), count);

        } else {
            neighbors =  getVertexNeighbors(v, count);
        }
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    const ui getSpatialIndex(VertexID id) const {
        return spatial_index_[id];
    }

    const ui getMaxSpatialIndex() const {
        return max_spatial_index_;
    }

    void buildCoreTable();

    void computeSpatialDistance();
};


#endif //SUBGRAPHMATCHING_GRAPH_H
