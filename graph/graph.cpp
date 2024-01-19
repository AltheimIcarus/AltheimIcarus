//
// Created by ssunah on 6/22/18.
//

#include "graph.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <map>
#include <utility/graphoperations.h>
#include <queue>

void Graph::BuildReverseIndex() {
    reverse_index_ = new ui[vertices_count_];
    reverse_index_offsets_= new ui[labels_count_ + 1];
    reverse_index_offsets_[0] = 0;

    ui total = 0;
    for (ui i = 0; i < labels_count_; ++i) {
        reverse_index_offsets_[i + 1] = total;
        total += labels_frequency_[i];
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        LabelID label = labels_[i];
        reverse_index_[reverse_index_offsets_[label + 1]++] = i;
    }
}

#if OPTIMIZED_LABELED_GRAPH == 1
void Graph::BuildNLF() {
    neighbors_by_labels_ = new VertexID[edges_count_ * 2];
    nlf_ = new std::unordered_map<LabelID, ui>[vertices_count_];
    nlo = new std::unordered_map<LabelID, std::pair<ui,ui>>[vertices_count_];
    
    // Label Pair Index [L(u)][L(v)] = [u']
    label_pairs.resize(labels_count_);

    for (ui i = 0; i < vertices_count_; ++i) {
        std::map<LabelID, std::vector<VertexID>> labelled_adjs;
        ui start = offsets_[i];
        ui end = offsets_[i+1];
        LabelID u_label = getVertexLabel(i);
        ui count;
        const VertexID * neighbors = getVertexNeighbors(i, count);

        for (ui j = 0; j < count; ++j) {
            VertexID v = neighbors[j];
            LabelID v_label = getVertexLabel(v);
            if (nlf_[i].find(v_label) == nlf_[i].end()) {
                nlf_[i][v_label] = 0;
            }

            nlf_[i][v_label] += 1;

            // Label Pair
            auto lp_iter = label_pairs[u_label].emplace(v_label, std::vector<ui>()); // create new label pair if not exist
            // if u' is not in set of LP[L(u')][L(v')]
            if (labelled_adjs.count(v_label)==0)
                lp_iter.first->second.push_back(i); // push u' into LP[L(u')][L(v')]

            // STORE: [label_u][label_v] = [u(s)], then [label_v][Label_u] = [v(s)] Equal to find left vertices that has adjs with label_v
            labelled_adjs[v_label].push_back(v);
        }

        ui ctr = start;
        for (auto it=labelled_adjs.begin(); it!=labelled_adjs.end(); ++it) {
            ui u_nlf = it->second.size();
            if (u_nlf == 0)
                continue;
            
            // nlo[u][label] = (start index, count) for list of adjs with same label id
            nlo[i].emplace(it->first, std::make_pair(ctr, u_nlf));
            
            // sort adjacent with label l to enable binary search
            std::sort(it->second.begin(), it->second.end());
            for (VertexID v : it->second) {
                neighbors_by_labels_[ctr++] = v;
            }

            if (ctr == end)
                break;
        }

        if (ctr != end) printf("NLO size %d\n", nlo[i].size());
    }
}

void Graph::BuildLabelOffset() {
    size_t labels_offset_size = (size_t)vertices_count_ * labels_count_ + 1;
    labels_offsets_ = new ui[labels_offset_size];
    std::fill(labels_offsets_, labels_offsets_ + labels_offset_size, 0);

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1],
            [this](const VertexID u, const VertexID v) -> bool {
                return labels_[u] == labels_[v] ? u < v : labels_[u] < labels_[v];
            });
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        LabelID previous_label = 0;
        LabelID current_label = 0;

        labels_offset_size = i * labels_count_;
        labels_offsets_[labels_offset_size] = offsets_[i];

        for (ui j = offsets_[i]; j < offsets_[i + 1]; ++j) {
            current_label = labels_[neighbors_[j]];

            if (current_label != previous_label) {
                for (ui k = previous_label + 1; k <= current_label; ++k) {
                    labels_offsets_[labels_offset_size + k] = j;
                }
                previous_label = current_label;
            }
        }

        for (ui l = current_label + 1; l <= labels_count_; ++l) {
            labels_offsets_[labels_offset_size + l] = offsets_[i + 1];
        }
    }
}

#endif

void Graph::loadGraphFromFile(const std::string &file_path) {
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    infile >> type >> vertices_count_ >> edges_count_;
    offsets_ = new ui[vertices_count_ +  1];
    offsets_[0] = 0;

    neighbors_ = new VertexID[edges_count_ * 2];
    labels_ = new LabelID[vertices_count_];
    labels_count_ = 0;
    max_degree_ = 0;

    LabelID max_label_id = 0;
    std::vector<ui> neighbors_offset(vertices_count_, 0);

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            VertexID id;
            LabelID  label;
            ui degree;
            infile >> id >> label >> degree;

            labels_[id] = label;
            offsets_[id + 1] = offsets_[id] + degree;

            if (degree > max_degree_) {
                max_degree_ = degree;
            }

            if (labels_frequency_.find(label) == labels_frequency_.end()) {
                labels_frequency_[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            labels_frequency_[label] += 1;
        }
        else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
            infile >> begin >> end;

            ui offset = offsets_[begin] + neighbors_offset[begin];
            neighbors_[offset] = end;

            offset = offsets_[end] + neighbors_offset[end];
            neighbors_[offset] = begin;

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }

    infile.close();
    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    BuildReverseIndex();
    //computeSpatialDistance();

#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        // BuildLabelOffset();
    }
#endif
}

void Graph::loadGraphFromVertexFile(const std::string &vertex_path, const std::string &edge_path) {
    std::ifstream infile(vertex_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the vertex file " << vertex_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    infile >> type >> vertices_count_ >> edges_count_;
    offsets_ = new ui[vertices_count_ +  1];
    offsets_[0] = 0;

    neighbors_ = new VertexID[edges_count_ * 2];
    labels_ = new LabelID[vertices_count_];
    labels_count_ = 0;
    max_degree_ = 0;

    LabelID max_label_id = 0;

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            VertexID id;
            LabelID  label;
            ui degree;
            infile >> id >> label >> degree;

            labels_[id] = label;
            offsets_[id + 1] = offsets_[id] + degree;

            if (degree > max_degree_) {
                max_degree_ = degree;
            }

            if (labels_frequency_.find(label) == labels_frequency_.end()) {
                labels_frequency_[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            labels_frequency_[label] += 1;
        }
    }

    infile.close();
    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    BuildReverseIndex();
    //computeSpatialDistance();

    loadGraphFromEdgeFile(edge_path);
}

void Graph::loadGraphFromEdgeFile(const std::string &file_path) {
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the edge file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    
    std::vector<ui> neighbors_offset(vertices_count_, 0);

    while (infile >> type) {
        if (type == 'e') { // Read edge only.
            VertexID begin;
            VertexID end;
            infile >> begin >> end;

            ui offset = offsets_[begin] + neighbors_offset[begin];
            neighbors_[offset] = end;

            offset = offsets_[end] + neighbors_offset[end];
            neighbors_[offset] = begin;

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }

    infile.close();
    
    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    //computeSpatialDistance();

#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        // BuildLabelOffset();
    }
#endif
}


void Graph::printGraphMetaData() {
    std::cout << "|V|: " << vertices_count_ << ", |E|: " << edges_count_ << ", |\u03A3|: " << labels_count_ << std::endl;
    std::cout << "Max Degree: " << max_degree_ << ", Max Label Frequency: " << max_label_frequency_ << std::endl;
}

void Graph::buildCoreTable() {
    core_table_ = new int[vertices_count_];
    GraphOperations::getKCore(this, core_table_);

    for (ui i = 0; i < vertices_count_; ++i) {
        if (core_table_[i] > 1) {
            core_length_ += 1;
        }
    }
}

void Graph::loadGraphFromFileCompressed(const std::string &degree_path, const std::string &edge_path,
                                        const std::string &label_path) {
    std::ifstream deg_file(degree_path, std::ios::binary);

    if (deg_file.is_open()) {
        std::cout << "Open degree file " << degree_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open degree file " << degree_path << " ." << std::endl;
        exit(-1);
    }

    auto start = std::chrono::high_resolution_clock::now();
    int int_size;
    deg_file.read(reinterpret_cast<char *>(&int_size), 4);
    deg_file.read(reinterpret_cast<char *>(&vertices_count_), 4);
    deg_file.read(reinterpret_cast<char *>(&edges_count_), 4);

    offsets_ = new ui[vertices_count_ + 1];
    ui* degrees = new unsigned int[vertices_count_];

    deg_file.read(reinterpret_cast<char *>(degrees), sizeof(int) * vertices_count_);


    deg_file.close();
    deg_file.clear();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Load degree file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    std::ifstream adj_file(edge_path, std::ios::binary);

    if (adj_file.is_open()) {
        std::cout << "Open edge file " << edge_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open edge file " << edge_path << " ." << std::endl;
        exit(-1);
    }

    start = std::chrono::high_resolution_clock::now();
    size_t neighbors_count = (size_t)edges_count_ * 2;
    neighbors_ = new ui[neighbors_count];

    offsets_[0] = 0;
    for (ui i = 1; i <= vertices_count_; ++i) {
        offsets_[i] = offsets_[i - 1] + degrees[i - 1];
    }

    max_degree_ = 0;

    for (ui i = 0; i < vertices_count_; ++i) {
        if (degrees[i] > 0) {
            if (degrees[i] > max_degree_)
                max_degree_ = degrees[i];
            adj_file.read(reinterpret_cast<char *>(neighbors_ + offsets_[i]), degrees[i] * sizeof(int));
            std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
        }
    }

    adj_file.close();
    adj_file.clear();

    delete[] degrees;

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Load adj file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;


    std::ifstream label_file(label_path, std::ios::binary);
    if (label_file.is_open())  {
        std::cout << "Open label file " << label_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open label file " << label_path << " ." << std::endl;
        exit(-1);
    }

    start = std::chrono::high_resolution_clock::now();

    labels_ = new ui[vertices_count_];
    label_file.read(reinterpret_cast<char *>(labels_), sizeof(int) * vertices_count_);

    label_file.close();
    label_file.clear();

    ui max_label_id = 0;
    for (ui i = 0; i < vertices_count_; ++i) {
        ui label = labels_[i];

        if (labels_frequency_.find(label) == labels_frequency_.end()) {
            labels_frequency_[label] = 0;
            if (label > max_label_id)
                max_label_id = label;
        }

        labels_frequency_[label] += 1;
    }

    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Load label file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    BuildReverseIndex();
    //computeSpatialDistance();
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Build reverse index file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        // BuildLabelOffset();
    }
#endif
}

void Graph::storeComparessedGraph(const std::string& degree_path, const std::string& edge_path,
                                  const std::string& label_path) {
    ui* degrees = new ui[vertices_count_];
    for (ui i = 0; i < vertices_count_; ++i) {
        degrees[i] = offsets_[i + 1] - offsets_[i];
    }

    std::ofstream deg_outputfile(degree_path, std::ios::binary);

    if (deg_outputfile.is_open()) {
        std::cout << "Open degree file " << degree_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot degree edge file " << degree_path << " ." << std::endl;
        exit(-1);
    }

    int int_size = sizeof(int);
    size_t vertex_array_bytes = ((size_t)vertices_count_) * 4;
    deg_outputfile.write(reinterpret_cast<const char *>(&int_size), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&vertices_count_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&edges_count_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(degrees), vertex_array_bytes);

    deg_outputfile.close();
    deg_outputfile.clear();

    delete[] degrees;

    std::ofstream edge_outputfile(edge_path, std::ios::binary);

    if (edge_outputfile.is_open()) {
        std::cout << "Open edge file " << edge_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot edge file " << edge_path << " ." << std::endl;
        exit(-1);
    }

    size_t edge_array_bytes = ((size_t)edges_count_ * 2) * 4;
    edge_outputfile.write(reinterpret_cast<const char *>(neighbors_), edge_array_bytes);

    edge_outputfile.close();
    edge_outputfile.clear();

    std::ofstream label_outputfile(label_path, std::ios::binary);

    if (label_outputfile.is_open()) {
        std::cout << "Open label file " << label_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot label file " << label_path << " ." << std::endl;
        exit(-1);
    }

    size_t label_array_bytes = ((size_t)vertices_count_) * 4;
    label_outputfile.write(reinterpret_cast<const char *>(labels_), label_array_bytes);

    label_outputfile.close();
    label_outputfile.clear();
}


void Graph::computeSpatialDistance() {
    spatial_index_ = new ui[vertices_count_];
    std::vector<bool> visited(vertices_count_, false);
    std::queue<ui> q;

    ui idx = 0;
    
    for (ui i=0; i<vertices_count_; ++i) {
        if (visited[i])
            continue;

        ui last_v = i;
        q.push(i);
        visited[i] = true;
        spatial_index_[i] = idx;
        // BFS Traversal
        while (!q.empty()) {
            VertexID u = q.front();
            q.pop();
            // all pair of connected vertices either has:
            //      difference of -1 (previous depth) or +1 (next depth),
            //      or same spatial index (same level of depth)

            ui deg;
            const ui* adjs = getVertexNeighbors(u, deg);
            
            for (ui j=0; j<deg; ++j) {
                ui v = adjs[j];
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                    spatial_index_[v] = spatial_index_[u] + 1;
                }
            }

            if (u == last_v) {
                idx += 1;
                // let the last vertex in next depth = last element in the queue
                last_v = q.back();
            }
        }

        // increment spatial distance by 2 to separate next disjoint graph
        idx += 2;
    }
    max_spatial_index_ = idx - 1;
}


void Graph::loadGraphFromFileBeta(const std::string &file_path) {
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    infile >> vertices_count_;
    offsets_ = new ui[vertices_count_ +  1];
    offsets_[0] = 0;
    labels_ = new LabelID[vertices_count_];
    
    std::vector<VertexID>tmp_edges;
    tmp_edges.reserve(vertices_count_ * (vertices_count_ - 1) / 2);
    std::string line;
    std::string token;

    VertexID u = 0;
    edges_count_ = 0;
    avg_deg_ = 0.0;

    while (u < vertices_count_ && std::getline(infile, line)) {
        labels_[u] = 0;

        std::stringstream ss(line);
        infile >> token; // degree

        int deg = std::stoi(token);
        avg_deg_ += deg;
        //printf("%d[%d]\t", u, deg);
        offsets_[u+1] = offsets_[u] + deg;
        while (deg > 0) {
            infile >> token;
            VertexID v = std::stoi(token);
            //printf("%d,", v);
            tmp_edges.push_back(v);
            if (v > u)
                edges_count_++;
            deg--;
        }
        u++;
        //printf("\n");
    }

    avg_deg_ /= (double) (vertices_count_ * 1.0);

    labels_count_ = 0;
    max_degree_ = 0;

    labels_count_ = 1;
    labels_frequency_[0] = vertices_count_;
    max_label_frequency_ = vertices_count_;

    neighbors_ = new VertexID[edges_count_ * 2];
    for (ui i = 0; i < vertices_count_; ++i) {
        ui deg = offsets_[i+1] - offsets_[i];
        if (deg > max_degree_) {
            max_degree_ = deg;
        }

        for (ui j = offsets_[i]; j < offsets_[i+1]; ++j) { // Read edge.
            neighbors_[j] = tmp_edges[j];
        }

        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    infile.close();

    BuildReverseIndex();
    //computeSpatialDistance();

#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        // BuildLabelOffset();
    }
#endif
}

