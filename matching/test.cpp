//
// Created by Roan Cheng Chi Qin on 21/03/23 DD/MM/YY.
//

#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <map>
#include <filesystem>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <iomanip>


#include "matchingcommand.h"
#include "graph/graph.h"
#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include "BuildTable.h"
#include "GenerateQueryPlan.h"
#include "EvaluateQuery.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

namespace fs = std::filesystem;

template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

template<typename Arg>
std::string doubleToString(Arg input, ui precision)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << input;
    return stream.str();
}

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}


#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))


bool test(const int algo, const Graph* data, const Graph* query,
    double& avg_candidate_count,
    double& avg_filtering_time_ns,
    double& avg_indexing_time_ns,
    double& avg_ordering_time_ns,
    double& avg_enumeration_time_ns,
    double& avg_query_time_ns,
    size_t& avg_memory_cost_bytes,
    size_t& avg_embedding_count,
    size_t& avg_call_count,
    double& total_time_ns)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    ui query_size = query->getVerticesCount();
    /* 
     * Filtering and Pruning
     */
    printf("filtering\n");

    ui** mapping = NULL;
    ui* candidate_count = NULL;
    bool successful_filter = false;
    ui* cfl_order = NULL;
    TreeNode* cfl_tree = NULL;
    ui* dpiso_order = NULL;
    TreeNode* dpiso_tree = NULL;
    TreeNode* ceci_tree = NULL;
    ui* ceci_order = NULL;
    TreeNode* l2Match_tree = NULL;
    ui* l2Match_order = NULL;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_candidate;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_candidate;


    start = std::chrono::high_resolution_clock::now();

    if (algo == 0)
    {
        // Novel l2Match
        successful_filter =  FilterVertices::l2MatchFilter(data, query, mapping, candidate_count, l2Match_order, l2Match_tree, TE_candidate, NTE_candidate);
    }
    else if (algo == 1)
    {
        // CECI 2019
        successful_filter =  FilterVertices::CECIFilter(data, query, mapping, candidate_count, ceci_order, ceci_tree, TE_candidate, NTE_candidate);
    }
    else if (algo == 2)
    {
        // CFL 2016
        successful_filter = FilterVertices::CFLFilter(data, query, mapping, candidate_count, cfl_order, cfl_tree);
    }
    else if (algo == 3)
    {
        // GQL 2008
        successful_filter = FilterVertices::GQLFilter(data, query, mapping, candidate_count);
    }
    else if (algo == 4)
    {
        // DP-ISO 2019
        successful_filter = FilterVertices::DPisoFilter(data, query, mapping, candidate_count, dpiso_order, dpiso_tree);
    }


    end = std::chrono::high_resolution_clock::now();
    double filtering_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    avg_filtering_time_ns += filtering_time_ns;
    
    if (!successful_filter)
    {
        std::cout << "No possible subgraph matched." << std::endl;
        return false;
    }

    ui total_candidates = 0;
    for (ui i=0; i<query_size; ++i)
    {
        total_candidates += candidate_count[i];
    }
    avg_embedding_count += total_candidates /(double) query_size;

    /* 
     * Indexing
     */
    printf("indexing\n");
    //std::cout << "BuildIndex()" << std::endl;

    Edges*** edge_mapping = NULL;
    double indexing_time_ns = 0.0;

    if (algo > 1)
    {
        edge_mapping = new Edges**[query->getVerticesCount()];
        for (ui i=0; i<query->getVerticesCount(); ++i)
        {
            edge_mapping[i] = new Edges*[query->getVerticesCount()];
        }

        start = std::chrono::high_resolution_clock::now();
        
        BuildTable::buildTables(data, query, mapping, candidate_count, edge_mapping);

        end = std::chrono::high_resolution_clock::now();
        indexing_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    avg_indexing_time_ns += indexing_time_ns;

    size_t memory_cost_bytes = 0;
    if (algo == 0)
    {
        // Novel l2Match_
        memory_cost_bytes = BuildTable::computeMemoryCostInBytes(query, candidate_count, l2Match_order, l2Match_tree, TE_candidate, NTE_candidate);
    }
    else if (algo == 1) {
        // CECI 2019
        memory_cost_bytes = BuildTable::computeMemoryCostInBytes(query, candidate_count, ceci_order, ceci_tree, TE_candidate, NTE_candidate);
    }
    else
    {
        // GQL 2008, CFL 2016, DP-ISO 2019
        memory_cost_bytes = BuildTable::computeMemoryCostInBytes(query, candidate_count, edge_mapping);
    }

    avg_memory_cost_bytes += memory_cost_bytes;
    

    /* 
     * Ordering
     */
    printf("ordering\n");
    //std::cout << "Ordering()" << std::endl;

    ui* matching_order = NULL;
    ui* pivots = NULL;
    ui** weight_array = NULL;
    
    start = std::chrono::high_resolution_clock::now();

    if (algo == 0)
    {
        // l2Match
        GenerateQueryPlan::generateCECIQueryPlan(query, l2Match_tree, l2Match_order, matching_order, pivots);
    }
    else if (algo == 1) {
        // CECI 2019
        GenerateQueryPlan::generateCECIQueryPlan(query, ceci_tree, ceci_order, matching_order, pivots);
    }
    else if (algo == 2)
    {
        // CFL 2016
        GenerateQueryPlan::generateCFLQueryPlan(data, query, edge_mapping, matching_order, pivots, cfl_tree, cfl_order, candidate_count);
    }
    else if (algo == 3)
    {
        // GQL 2008
        GenerateQueryPlan::generateGQLQueryPlan(data, query, candidate_count, matching_order, pivots);
    }
    else if (algo == 4)
    {
        // DP-ISO 2019
        GenerateQueryPlan::generateDSPisoQueryPlan(query, edge_mapping, matching_order, pivots, dpiso_tree, dpiso_order, candidate_count, weight_array);
    }

    end = std::chrono::high_resolution_clock::now();
    double ordering_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    avg_ordering_time_ns += ordering_time_ns;

    /* 
     * Enumeration
     */
    printf("enumeration\n");
    //std::cout << "Enumerate()" << std::endl;

    size_t call_count = 0;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t embedding_count = 0;
    TimeOutException* timeout = new TimeOutException();
    
    start = std::chrono::high_resolution_clock::now();

    if (algo == 0)
    {
        // l2Match
        embedding_count = EvaluateQuery::exploreCECIStyle(data, query, l2Match_tree, mapping, candidate_count, TE_candidate, NTE_candidate, l2Match_order, output_limit, call_count, timeout);
    }
    else if (algo == 1) {
        // CECI 2019
        embedding_count = EvaluateQuery::exploreCECIStyle(data, query, ceci_tree, mapping, candidate_count, TE_candidate, NTE_candidate, ceci_order, output_limit, call_count, timeout);
    }
    else if (algo == 2)
    {
        // CFL 2016
        embedding_count = EvaluateQuery::exploreGraph(data, query, edge_mapping, mapping, candidate_count, matching_order, pivots, output_limit, call_count, timeout);
    }
    else if (algo == 3)
    {
        // GQL 2008
        embedding_count = EvaluateQuery::exploreGraphQLStyle(data, query, mapping, candidate_count, matching_order, output_limit, call_count, timeout);
    }
    else if (algo == 4)
    {
        // DP-ISO 2019
        embedding_count = EvaluateQuery::exploreDPisoStyle(data, query, dpiso_tree, edge_mapping, mapping, candidate_count, weight_array, dpiso_order, output_limit, call_count, timeout);
    }

    end = std::chrono::high_resolution_clock::now();
    double enumeration_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    avg_enumeration_time_ns += enumeration_time_ns;
    total_time_ns += filtering_time_ns + indexing_time_ns + ordering_time_ns + enumeration_time_ns;
    avg_query_time_ns += total_time_ns;
    avg_embedding_count += embedding_count;
    avg_call_count += call_count;

    /* 
     * Release memory
     */
    delete[] cfl_order;
    delete[] cfl_tree;
    delete[] dpiso_order;
    delete[] dpiso_tree;
    delete[] ceci_order;
    delete[] ceci_tree;
    delete[] l2Match_order;
    delete[] l2Match_tree;
    delete[] matching_order;
    delete[] pivots;
    for (ui i=0; i<query->getVerticesCount(); ++i)
    {
        delete[] mapping[i];
    }
    delete[] mapping;
    delete[] candidate_count;
    
    if (edge_mapping != NULL)
    {
        for (ui i=0; i<query->getVerticesCount(); ++i)
        {
            for (ui j=0; j<query->getVerticesCount(); ++j)
                delete edge_mapping[i][j];
            
            delete[] edge_mapping[i];
        }
        delete[] edge_mapping;
    }
    if (weight_array != NULL)
    {
        for (ui i=0; i<query->getVerticesCount(); ++i)
            delete[] weight_array[i];

        delete[] weight_array;
    }

    return true;
}

bool valid_test_case(std::string dataset, int test_size)
{
    std::string human_dataset = "human";
    if (dataset == human_dataset)
    {
        if (test_size==24 || test_size==32)
            return false;
    }
    else
    {
        if (test_size==12 || test_size==20)
            return false;
    }
    return true;
}

std::map<int, std::string> getTestFiles(int query_vertices_count, std::string test_density, std::string test_dataset, std::string query_directory)
{
    std::string test_size = std::to_string(query_vertices_count); // |V(Q)|
    
    std::map<int, std::string> query_files;
    std::string query_file_prefix = test_density + test_size + "_";
    auto query_file_prefix_cstr = query_file_prefix.c_str();

    if (!valid_test_case(test_dataset, query_vertices_count))
        return query_files;

    // walk the query graphs directory
    for (auto dir = fs::recursive_directory_iterator(query_directory);
        dir != fs::recursive_directory_iterator();
        ++dir)
    {
        if (dir->is_regular_file() && dir->path().extension().string()==".graph")
        {
            std::string filename = dir->path().filename().string();
            int test_id = stoi(split(split(filename, '_').back(), '.')[0]); // get query graph ID
            // check for matching filename prefix
            if (strncmp(filename.c_str(), query_file_prefix_cstr, strlen(query_file_prefix_cstr))==0 
                && test_id < 201 )
            {
                query_files[test_id] = (dir->path().string());
            }
        }
    }

    // Return test case files
    return query_files;
}


int main_test() {
    const int ALGO_CNT = 5;
    // Import dataset
    /* 
     * Import Graph
     */

    // DATASETS
    std::vector<std::string> data_files{
        "../../test/dataset/dblp/data_graph/dblp.graph",
        "../../test/dataset/human/data_graph/human.graph",
        "../../test/dataset/patents/data_graph/patents.graph",
        "../../test/dataset/youtube/data_graph/youtube.graph",
        "../../test/dataset/yeast/data_graph/yeast.graph"
    };
    std::vector<std::string> query_dirs{
        "../../test/dataset/dblp/query_graph/",
        "../../test/dataset/human/query_graph/",
        "../../test/dataset/patents/query_graph/",
        "../../test/dataset/youtube/query_graph/",
        "../../test/dataset/yeast/query_graph/"
    };
    std::vector<std::pair<int, std::string>> test_sizes{
        {4, "query_dense_"},
        {8, "query_dense_"},
        {8, "query_sparse_"},
        {12, "query_dense_"},    // human dataset only
        {12, "query_sparse_"},   // human dataset only
        {16, "query_dense_"},
        {16, "query_sparse_"},
        {20, "query_dense_"},    // human dataset only
        {20, "query_sparse_"},   // human dataset only
        {24, "query_dense_"},    // except human dataset
        {24, "query_sparse_"},   // except human dataset
        {32, "query_dense_"},    // except human dataset
        {32, "query_sparse_"}   // except human dataset
    };

    {
        int dataset = 0;
        int types = 5;
        int n = 101;
        std::string data_file(data_files[dataset]);
        Graph* data = new Graph(true);
        data->loadGraphFromFile(data_file);

        std::string test_size = std::to_string(test_sizes[types].first); // |V(Q)|
        std::string test_density(test_sizes[types].second); // density
        std::string test_dataset = split(data_file,'/')[2]; //dataset name, eg: dblp
        std::string query_file_prefix = std::string(query_dirs[dataset]) + test_density + test_size + "_" + std::to_string(n) + ".graph";
        auto query_file_prefix_cstr = query_file_prefix.c_str();
        Graph* query = new Graph(true);
        query->loadGraphFromFile(query_file_prefix);
        query->buildCoreTable();

        std::vector<std::string> algos = {
            "1) LPF;",
            "2) CECI;",
            "3) CFL;",
            "4) GQL;",
            "5) DP-ISO;"
        };

        std::string result;
        result = (std::string("Algorithm;")
            + "AVG Candidate Count;" + "AVG Filtering Time (s);"
            + "AVG Indexing Time (s);" + "AVG Ordering Time (s);"
            + "AVG Enumeration Time (s);" + "AVG Query Time (s);"
            + "AVG Memory Cost (MB);" + "AVG Embedding Count;"
            + "AVG Call Count;\n"
        );
        std::cout << result;
        
        // loop through each algorithms
        for (int a_idx=1; a_idx<ALGO_CNT; ++a_idx)
        {
            double candidate_count = 0.0;
            double filtering_time_s = 0.0;
            double indexing_time_s = 0.0;
            double ordering_time_s = 0.0;
            double enum_time_s = 0.0;
            double query_time_s = 0.0;
            size_t mem_cost_MB = 0;
            size_t embedding_count = 0;
            size_t call_count = 0;
            double total_query_time_s = 0.0;

            bool valid_result = test(
                a_idx, data, query,
                candidate_count,
                filtering_time_s,
                indexing_time_s,
                ordering_time_s,
                enum_time_s,
                query_time_s,
                mem_cost_MB,
                embedding_count,
                call_count,
                total_query_time_s
            );

            result = (algos[a_idx] + ";"
                + doubleToString(candidate_count, 8) + ";" + doubleToString(filtering_time_s/(double)1000000000, 6) + ";"
                + doubleToString(indexing_time_s/(double)1000000000, 6) + ";" + doubleToString(ordering_time_s/(double)1000000000, 6) + ";"
                + doubleToString(enum_time_s/(double)1000000000, 6) + ";" + doubleToString(query_time_s/(double)1000000000, 6) + ";"
                + doubleToString(mem_cost_MB/(double)(1024 * 1024), 6) + ";" + doubleToString(embedding_count, 8) + ";"
                + doubleToString(call_count, 6) + ";\n"
            );
            std::cout << result;
        }

        delete data;
        delete query;
        return 0;
    }

    for (int i=0; i<data_files.size(); ++i) {
        // SET working directory
        std::string data_file(data_files[i]);
        std::string query_dir(query_dirs[i]);

        // Import data graph
        auto start = std::chrono::high_resolution_clock::now();

        Graph* data = new Graph(true);
        data->loadGraphFromFile(data_file);

        auto end = std::chrono::high_resolution_clock::now();
        double import_graph_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        data->printGraphMetaData();

        // for each combination of (|V(Q)|,DENSITY)
        for (auto test_case : test_sizes) {
            int test_size_int = test_case.first; // |V(Q)|
            std::string test_size = std::to_string(test_size_int); // |V(Q)|
            std::string test_density(test_case.second); // density
            std::string test_dataset = split(data_file,'/')[2]; //dataset name, eg: dblp

            std::map<int, std::string> query_files = getTestFiles(test_size_int, test_density, test_dataset, query_dir);
            if (!query_files.empty())
                continue;

            // Import test inputs
            int test_count = query_files.size();

            // Competing algorithm names
            std::vector<std::string> algos = {
                "1) LDF;",
                "2) LPF;",
                "3) CECI;",
                "4) CFL;",
                "5) GQL;",
                "6) DP-ISO;"
            };

            // results placeholder
            double* avg_candidate_count = new double[ALGO_CNT];
            double* avg_filtering_time_ns = new double[ALGO_CNT];

            double* avg_indexing_time_ns = new double[ALGO_CNT];
            
            double* avg_ordering_time_ns = new double[ALGO_CNT];
            
            double* avg_enumeration_time_ns = new double[ALGO_CNT];
            
            double* avg_query_time_ns = new double[ALGO_CNT];
            size_t* avg_memory_cost_bytes = new size_t[ALGO_CNT];
            size_t* avg_embedding_count = new size_t[ALGO_CNT];
            size_t* avg_call_count = new size_t[ALGO_CNT];
            double* total_query_time_ns = new double[ALGO_CNT];

            // loop through each query graph
            for(auto query_file : query_files) {
                Graph* query = new Graph(true);
                query->loadGraphFromFile(query_file.second);
                query->buildCoreTable();
                //query->printGraphMetaData();

                // loop through each algorithms
                for (int a_idx=0; a_idx<ALGO_CNT; ++a_idx) {
                    bool valid_result = test(
                        a_idx, data, query,
                        avg_candidate_count[a_idx],
                        avg_filtering_time_ns[a_idx],
                        avg_indexing_time_ns[a_idx],
                        avg_ordering_time_ns[a_idx],
                        avg_enumeration_time_ns[a_idx],
                        avg_query_time_ns[a_idx],
                        avg_memory_cost_bytes[a_idx],
                        avg_embedding_count[a_idx],
                        avg_call_count[a_idx],
                        total_query_time_ns[a_idx]
                    );

                    if (!valid_result)
                        continue;
                }
            }

            // Write results to a CSV file
            vector<std::string> results;
            results.emplace_back(std::string("Algorithm;")
                + "AVG Candidate Count;" + "AVG Filtering Time (s);"
                + "AVG Indexing Time (s);" + "AVG Ordering Time (s);"
                + "AVG Enumeration Time (s);" + "AVG Query Time (s);"
                + "AVG Memory Cost (MB);" + "AVG Embedding Count;"
                + "AVG Call Count;\n"
            );

            std::string root_dir = "result";
            fs::create_directory(root_dir);
            std::string output_file = root_dir + "/" + test_dataset + "_" + test_size + "_" + test_density + ".csv";

            
            // loop through each algorithms
            for (int a_idx=0; a_idx<ALGO_CNT; ++a_idx)
            {
                double candidate_count = avg_candidate_count[a_idx] / (double)test_count;
                double filtering_time_s = avg_filtering_time_ns[a_idx] / (double)test_count /(double)1000000000;
                double indexing_time_s = avg_indexing_time_ns[a_idx] / (double)test_count /(double)1000000000;
                double ordering_time_s = avg_ordering_time_ns[a_idx] / (double)test_count /(double)1000000000;
                double enum_time_s = avg_enumeration_time_ns[a_idx] / (double)test_count /(double)1000000000;
                double query_time_s = avg_query_time_ns[a_idx] / (double)test_count /(double)1000000000;
                double mem_cost_MB = avg_memory_cost_bytes[a_idx] / test_count /(double)(1024 * 1024);
                double embedding_count = avg_embedding_count[a_idx] / (double)test_count;
                double call_count = avg_call_count[a_idx] / (double)test_count;

                results.emplace_back(algos[a_idx] + ";"
                    + doubleToString(candidate_count, 8) + ";" + doubleToString(filtering_time_s, 6) + ";"
                    + doubleToString(indexing_time_s, 6) + ";" + doubleToString(ordering_time_s, 6) + ";"
                    + doubleToString(enum_time_s, 6) + ";" + doubleToString(query_time_s, 6) + ";"
                    + doubleToString(mem_cost_MB, 6) + ";" + doubleToString(embedding_count, 8) + ";"
                    + doubleToString(call_count, 6) + ";\n"
                );
            }

            std::ofstream s_out(output_file, std::fstream::out | std::fstream::trunc);
            for (auto buffer : results)
            {
                s_out << buffer << std::endl;
            }
            s_out.close();
            std::cout << "result details saved to: " << output_file << PRINT_SEPARATOR;
        }
    }

    return 0;
}


int main(int argc, char** argv) {
    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string input_filter_type = command.getFilterType();
    std::string dataset = split(input_query_graph_file, '/')[4];
    bool is_debug = true;
    if (argc == 8) {
        std::cout << argv[7] << std::endl;
        is_debug = true;
    }
    /**
     * Output the command line information.
     */
    std::cout << "Command Line:" << std::endl;
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << std::endl;
    std::cout << "\tFilter Type: " << input_filter_type << std::endl;
    std::cout << "\tDataset: " << dataset << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Load input graphs.
     * ./test.out -d ../../test/dataset/phase/phase-30-150-0.34-0.36-0.42-1-target -q ../../test/dataset/phase/phase-30-150-0.34-0.36-0.42-1-pattern -filter CFL
     */
    //std::cout << "Load graphs..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);
    //query_graph->loadGraphFromFileBeta(input_query_graph_file);
    query_graph->buildCoreTable();

    Graph* data_graph = new Graph(true);

    data_graph->loadGraphFromFile(input_data_graph_file);
    
    // label dataset
    // std::string edge_file = "../../test/dataset/label/data_graph/10000_0.05_20.edges";
    // data_graph->loadGraphFromVertexFile(input_data_graph_file, edge_file);
    
    // ERP dataset
    // data_graph->loadGraphFromFileBeta(input_data_graph_file);

    auto end = std::chrono::high_resolution_clock::now();

    //double load_graphs_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    int query_vertices_count = query_graph->getVerticesCount();
    int query_edges_count = query_graph->getEdgesCount();
    double query_density = (2.0 * query_graph->getEdgesCount()) / (double)(query_graph->getVerticesCount() * (query_graph->getVerticesCount() - 1.0));
    double data_density = (2.0 * data_graph->getEdgesCount()) / (double)(data_graph->getVerticesCount() * (data_graph->getVerticesCount() - 1.0));
    //double data_average_deg = data_graph->getAverageDegree();
    //std::cout << "-----" << std::endl;
    //std::cout << "Query Graph Meta Information" << std::endl;
    //query_graph->printGraphMetaData();
    //std::cout << "-----" << std::endl;
    //data_graph->printGraphMetaData();

    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Start queries.
     */

    //std::cout << "Start queries..." << std::endl;
    //std::cout << "-----" << std::endl;
    std::cout << "Filter candidates... ";

    start = std::chrono::high_resolution_clock::now();

    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* tso_order = NULL;
    TreeNode* tso_tree = NULL;
    ui* cfl_order = NULL;
    TreeNode* cfl_tree = NULL;
    ui* dpiso_order = NULL;
    TreeNode* dpiso_tree = NULL;
    TreeNode* ceci_tree = NULL;
    ui* ceci_order = NULL;
    TreeNode* l2Match_tree = NULL;
    ui* l2Match_order = NULL;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    bool build_edge_index = true;
    bool filter_pass = false;
    
    if (input_filter_type == "l2MatchJR") { // complete version
        build_edge_index = false;
        filter_pass = FilterVertices::l2MatchFilter(data_graph, query_graph, candidates, candidates_count, l2Match_order, l2Match_tree, TE_Candidates, NTE_Candidates);
        //std::cout << l2Match_order[0] << " [" << candidates_count[l2Match_order[0]] << ", " << candidates[l2Match_order[0]][0] << "]\n";
        //FilterVertices::l2MatchFilterB(data_graph, query_graph, candidates, candidates_count, l2Match_order, l2Match_tree);
    } else if (input_filter_type == "l2MatchA") { // w/o LPF
        build_edge_index = false;
        FilterVertices::l2MatchFilterA(data_graph, query_graph, candidates, candidates_count, l2Match_order, l2Match_tree, TE_Candidates, NTE_Candidates);
    } else if (input_filter_type == "l2MatchB") { // w/o BCPRefine
        build_edge_index = false;
        FilterVertices::l2MatchFilterB(data_graph, query_graph, candidates, candidates_count, l2Match_order, l2Match_tree, TE_Candidates, NTE_Candidates);
    }  else if (input_filter_type == "l2Match") { // w/o JR
        build_edge_index = true;
        FilterVertices::l2MatchFilterNew(data_graph, query_graph, candidates, candidates_count, l2Match_order, l2Match_tree);
    } else if (input_filter_type == "GQL") {
        build_edge_index = true;
        filter_pass = FilterVertices::GQLFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "CFL") {
        build_edge_index = true;
        filter_pass = FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, cfl_order, cfl_tree);
    } else if (input_filter_type == "DPiso") {
        build_edge_index = true;
        filter_pass = FilterVertices::DPisoFilter(data_graph, query_graph, candidates, candidates_count, dpiso_order, dpiso_tree);
    } else if (input_filter_type == "CECI") {
        build_edge_index = false;
        filter_pass = FilterVertices::CECIFilter(data_graph, query_graph, candidates, candidates_count, ceci_order, ceci_tree, TE_Candidates, NTE_Candidates);
    } else if (input_filter_type == "TEST") {
        build_edge_index = true;
        FilterVertices::l2MatchFilterNew(data_graph, query_graph, candidates, candidates_count,
                                            l2Match_order, l2Match_tree);
    } else {
        std::cout << "The specified filter type '" << input_filter_type << "' is not supported." << std::endl;
        exit(-1);
    }

    // Sort the candidates to support the set intersections
    if (build_edge_index)
        FilterVertices::sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    end = std::chrono::high_resolution_clock::now();
    double filter_vertices_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << NANOSECTOSEC(filter_vertices_time_in_ns) << std::endl;

    size_t total_candidates_count = 0;

    for(auto i=0; i<query_graph->getVerticesCount(); ++i)
        total_candidates_count += candidates_count[i];
    std::cout << "total_candidates_count: " << total_candidates_count << std::endl;
    
    //std::cout << "-----" << std::endl;
    std::cout << "Build indices... ";

    start = std::chrono::high_resolution_clock::now();

    Edges ***edge_matrix = NULL;
    if (build_edge_index) {
        edge_matrix = new Edges **[query_graph->getVerticesCount()];
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
        }

        BuildTable::buildTables(data_graph, query_graph, candidates, candidates_count, edge_matrix);
    }

    end = std::chrono::high_resolution_clock::now();
    double build_table_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << NANOSECTOSEC(build_table_time_in_ns) << std::endl;

    std::cout << "Compute indices memory cost... ";

    size_t memory_cost_in_bytes = 0;
    if (input_filter_type == "l2MatchJR") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, l2Match_order, l2Match_tree,
                TE_Candidates, NTE_Candidates);
        //memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, edge_matrix);
        //BuildTable::printTableCardinality(query_graph, l2Match_tree, l2Match_order, TE_Candidates, NTE_Candidates);
    } else if (input_filter_type == "l2MatchA") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, l2Match_order, l2Match_tree,
                TE_Candidates, NTE_Candidates);
        //BuildTable::printTableCardinality(query_graph, l2Match_tree, l2Match_order, TE_Candidates, NTE_Candidates);
    } else if (input_filter_type == "l2MatchB") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, l2Match_order, l2Match_tree,
                TE_Candidates, NTE_Candidates);
        //BuildTable::printTableCardinality(query_graph, l2Match_tree, l2Match_order, TE_Candidates, NTE_Candidates);
    } else if (input_filter_type == "l2Match") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, edge_matrix);
        //BuildTable::printTableCardinality(query_graph, l2Match_tree, l2Match_order, TE_Candidates, NTE_Candidates);
    } else if (input_filter_type == "CECI") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, ceci_order, ceci_tree,
                TE_Candidates, NTE_Candidates);
        //BuildTable::printTableCardinality(query_graph, ceci_tree, ceci_order, TE_Candidates, NTE_Candidates);
    } else {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, edge_matrix);
        //BuildTable::printTableCardinality(query_graph, edge_matrix);
    }
    std::cout << BYTESTOMB(memory_cost_in_bytes) << std::endl;
    //std::cout << "-----" << std::endl;
    std::cout << "Generate a matching order... ";

    start = std::chrono::high_resolution_clock::now();

    ui* matching_order = NULL;
    ui* pivots = NULL;
    ui** weight_array = NULL;

    if (input_filter_type == "l2MatchJR") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, l2Match_tree, l2Match_order, matching_order, pivots);
        //GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, l2Match_tree, l2Match_order, candidates_count);
    } else if (input_filter_type == "l2MatchA") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, l2Match_tree, l2Match_order, matching_order, pivots);
        //GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, l2Match_tree, l2Match_order, candidates_count);
    } else if (input_filter_type == "l2MatchB") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, l2Match_tree, l2Match_order, matching_order, pivots);
        //GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, l2Match_tree, l2Match_order, candidates_count);
    } else if (input_filter_type == "l2Match") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, l2Match_tree, l2Match_order, matching_order, pivots);
        //GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, l2Match_tree, l2Match_order, candidates_count);
    } else if (input_filter_type == "GQL") {
        GenerateQueryPlan::generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);
    } else if (input_filter_type == "CFL") {
        if (cfl_tree == NULL) {
            int level_count;
            ui* level_offset;
            GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, cfl_tree, cfl_order, level_count, level_offset);
            delete[] level_offset;
        }
        GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, cfl_tree, cfl_order, candidates_count);
    } else if (input_filter_type == "DPiso") {
        if (dpiso_tree == NULL) {
            GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, dpiso_tree, dpiso_order);
        }

        GenerateQueryPlan::generateDSPisoQueryPlan(query_graph, edge_matrix, matching_order, pivots, dpiso_tree, dpiso_order,
                                                    candidates_count, weight_array);
    } else if (input_filter_type == "CECI") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, ceci_tree, ceci_order, matching_order, pivots);
    } else if (input_filter_type == "TEST") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, cfl_tree, cfl_order, matching_order, pivots);
        //GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, l2Match_tree, l2Match_order, candidates_count);
    }
    else {
        std::cout << "The specified order type '" << input_filter_type << "' is not supported." << std::endl;
    }

    end = std::chrono::high_resolution_clock::now();
    double generate_query_plan_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << NANOSECTOSEC(generate_query_plan_time_in_ns) << std::endl;

    //GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, matching_order, pivots);
    //GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);
    
    //std::cout << "-----" << std::endl;
    std::cout << "Enumerate... ";
    size_t output_limit = std::numeric_limits<size_t>::max();
    output_limit = 1000000;
    auto tokens = split(input_query_graph_file,'/');
    ui query_size = query_graph->getVerticesCount();
    //if (query_size > 8 || (tokens[4]=="human" && query_size > 4))
    //    output_limit = 1000000;
    size_t embedding_count = 0;
    size_t call_count = 0;
    double enumeration_time_in_ns = 0.0;


#if ENABLE_QFLITER == 1
    EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

    //goto ENDE;
    try {
        if (!filter_pass) throw filter_pass;
        //const auto time_limit = start + std::chrono::duration<double, std::nano>(300000000000); // end time = current time + 5 minutes
        TimeOutException* timeout_e = new TimeOutException();
        std::thread([&timeout_e]{
            std::this_thread::sleep_for(300s); // 5 minutes
            timeout_e->is_throwable = true;
        }).detach();

        start = std::chrono::high_resolution_clock::now();

        if (input_filter_type == "l2MatchJR") {
            embedding_count = EvaluateQuery::exploreL2Match(data_graph, query_graph, l2Match_tree, candidates, candidates_count, TE_Candidates,
                    NTE_Candidates, l2Match_order, output_limit, call_count, timeout_e);
            //embedding_count = EvaluateQuery::exploreGraph(data_graph, query_graph, edge_matrix, candidates,
            //                                              candidates_count, matching_order, pivots, output_limit, call_count);
        } else if (input_filter_type == "l2MatchA") {
            embedding_count = EvaluateQuery::exploreL2Match(data_graph, query_graph, l2Match_tree, candidates, candidates_count, TE_Candidates,
                    NTE_Candidates, l2Match_order, output_limit, call_count, timeout_e);
        } else if (input_filter_type == "l2MatchB") {
            embedding_count = EvaluateQuery::exploreL2Match(data_graph, query_graph, l2Match_tree, candidates, candidates_count, TE_Candidates,
                    NTE_Candidates, l2Match_order, output_limit, call_count, timeout_e);
        } else if (input_filter_type == "l2Match") { // w/o JR
            //embedding_count = EvaluateQuery::exploreL2MatchC(data_graph, query_graph, l2Match_tree, candidates, candidates_count, TE_Candidates,
            //        NTE_Candidates, l2Match_order, output_limit, call_count, timeout_e);
            embedding_count = EvaluateQuery::exploreL2MatchNew(data_graph, query_graph, edge_matrix, l2Match_tree, candidates, candidates_count,
                    l2Match_order, pivots, output_limit, call_count, timeout_e);
        } else if (input_filter_type == "CFL") {
            embedding_count = EvaluateQuery::exploreGraph(data_graph, query_graph, edge_matrix, candidates,
                                                        candidates_count, matching_order, pivots, output_limit, call_count, timeout_e);
        } else if (input_filter_type == "GQL") {
            embedding_count = EvaluateQuery::exploreGraphQLStyle(data_graph, query_graph, candidates, candidates_count,
                                                                matching_order, output_limit, call_count, timeout_e);
        }
        else if (input_filter_type == "DPiso") {
            embedding_count = EvaluateQuery::exploreDPisoStyle(data_graph, query_graph, dpiso_tree,
                                                            edge_matrix, candidates, candidates_count,
                                                            weight_array, dpiso_order, output_limit,
                                                            call_count, timeout_e);
        }
        else if (input_filter_type == "CECI") {
            embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph, ceci_tree, candidates, candidates_count, TE_Candidates,
                    NTE_Candidates, ceci_order, output_limit, call_count, timeout_e);
        }
        else if (input_filter_type == "TEST") {
            embedding_count = EvaluateQuery::exploreGraphJR(data_graph, query_graph, edge_matrix, candidates,
                                                        candidates_count, cfl_order, pivots, output_limit, call_count, timeout_e);
        }
        else {
            std::cout << "The specified engine type '" << input_filter_type << "' is not supported." << std::endl;
            exit(-1);
        }

        end = std::chrono::high_resolution_clock::now();
        enumeration_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        if (timeout_e->is_throwable)
            throw (*timeout_e);
    
    } catch (bool e) {
        enumeration_time_in_ns = 0.0;
    } catch (...) {
        enumeration_time_in_ns = 300000000000.0;
        std::cout << NANOSECTOSEC(enumeration_time_in_ns) << " : " << embedding_count << std::endl;
        std::cerr << "Timeout after 5 minutes" << std::endl;
    }
    std::cout << NANOSECTOSEC(enumeration_time_in_ns) << " : " << embedding_count << std::endl;

    //ENDE:
    
#ifdef DISTRIBUTION
    std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
    outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
    delete[] EvaluateQuery::distribution_count_;
#endif

    //std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
    /**
     * Release the allocated memories.
     */
    delete[] candidates_count;
    delete[] tso_order;
    delete[] tso_tree;
    delete[] cfl_order;
    delete[] cfl_tree;
    delete[] dpiso_order;
    delete[] dpiso_tree;
    delete[] ceci_order;
    delete[] ceci_tree;
    delete[] l2Match_order;
    delete[] l2Match_tree;
    delete[] matching_order;
    delete[] pivots;
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        delete[] candidates[i];
    }
    delete[] candidates;

    if (edge_matrix != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }
    if (weight_array != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            delete[] weight_array[i];
        }
        delete[] weight_array;
    }

    delete query_graph;
    delete data_graph;

    /**
     * End.
     */
    std::cout << "--------------------------------------------------------------------" << std::endl;
    double preprocessing_time_in_ns = filter_vertices_time_in_ns + build_table_time_in_ns + generate_query_plan_time_in_ns;
    double total_time_in_ns = preprocessing_time_in_ns + enumeration_time_in_ns;

    // Write results to a CSV file
    std::string root_dir = "result/phase";
    fs::create_directory(root_dir);
    
    std::string output_file = "";
    std::vector<std::string> results;

    if (is_debug) {
        double data_average_deg = data_graph->getAverageDegree();
        output_file = root_dir + "/result.csv";
        if (!fs::exists(output_file))
        {
            results.emplace_back(std::string("No.;") + "Algorithm;"
                + "Candidate Count;" + "Filtering Time (s);"
                + "Indexing Time (s);" + "Ordering Time (s);"
                + "Enumeration Time (s);" + "Query Time (s);"
                + "Embedding Count;"
                + "Call Count;" + "|V|;" + "|E|;" + "d;" + "|L(D)|;" + "p(Q);" + "p(D)"
            );
        }

        std::string output_filename = split(input_query_graph_file, '/')[6];
        auto test_parameters = split(output_filename, '_');
        std::string label_size = split(test_parameters[2], '.')[0];
        std::string query_p = test_parameters[1];

        std::string data_p = split(input_data_graph_file, '/')[6];
        data_p = split(data_p, '_')[1];

        results.emplace_back(output_filename + ";" + input_filter_type + ";"
            + std::to_string(total_candidates_count) + ";" + doubleToString(NANOSECTOSEC(filter_vertices_time_in_ns), 6) + ";"
            + doubleToString(NANOSECTOSEC(build_table_time_in_ns), 6) + ";" + doubleToString(NANOSECTOSEC(generate_query_plan_time_in_ns), 6) + ";"
            + doubleToString(NANOSECTOSEC(enumeration_time_in_ns), 6) + ";" + doubleToString(NANOSECTOSEC(total_time_in_ns), 6) + ";"
            + std::to_string(embedding_count) + ";"
            + std::to_string(call_count) + ";"
            + std::to_string(query_vertices_count) + ";"
            + std::to_string(query_edges_count) + ";"
            + doubleToString(query_density, 6) + ";"
            + label_size + ";" + query_p + ";" + data_p
        );
    } else {
        output_file = root_dir + "/LPF.csv";
        if (!fs::exists(output_file))
        {
            // append header if file not exists
            results.emplace_back(std::string("No.;") + "Algorithm;"
                + "Candidate Count;" + "Filtering Time (s);"
                + "Indexing Time (s);" + "Ordering Time (s);"
                + "Enumeration Time (s);" + "Query Time (s);"
                + "Embedding Count;"
                + "Call Count;" + "|V|;" + "|E|;" + "d;" + "dataset"
            );
        }

        //auto tokens = split(input_query_graph_file,'/');
        auto test_tokens = split(tokens[6], '_');
        std::string test_dataset = tokens[4]; //dataset name, eg: dblp
        std::string test_size(test_tokens[2]); // |V(Q)|
        std::string test_density(test_tokens[1]); // density: dense or sparse
        std::string test_index(split(test_tokens[3], '.')[0]); // 1-200 
        //std::string output_file = root_dir + "/" + test_dataset + ".csv";

        results.emplace_back(test_size+test_density+test_index + ";" + input_filter_type + ";"
            + std::to_string(total_candidates_count) + ";" + doubleToString(NANOSECTOSEC(filter_vertices_time_in_ns), 6) + ";"
            + doubleToString(NANOSECTOSEC(build_table_time_in_ns), 6) + ";" + doubleToString(NANOSECTOSEC(generate_query_plan_time_in_ns), 6) + ";"
            + doubleToString(NANOSECTOSEC(enumeration_time_in_ns), 6) + ";" + doubleToString(NANOSECTOSEC(total_time_in_ns), 6) + ";"
            + std::to_string(embedding_count) + ";"
            + std::to_string(call_count) + ";"
            + std::to_string(query_vertices_count) + ";"
            + std::to_string(query_edges_count) + ";"
            + doubleToString(query_density, 6) + ";"
            + test_dataset
        );

    }

    

    std::ofstream s_out(output_file, std::fstream::out | std::fstream::app);
    for (auto row : results)
        s_out << row << std::endl;
    
    s_out.close();
    //std::cout << "result details saved to: " << output_file << PRINT_SEPARATOR;

    // printf("Filter vertices time (seconds): %.4lf\n", NANOSECTOSEC(filter_vertices_time_in_ns));
    // printf("Build table time (seconds): %.4lf\n", NANOSECTOSEC(build_table_time_in_ns));
    // printf("Generate query plan time (seconds): %.4lf\n", NANOSECTOSEC(generate_query_plan_time_in_ns));
    // printf("Enumerate time (seconds): %.4lf\n", NANOSECTOSEC(enumeration_time_in_ns));
    // printf("Preprocessing time (seconds): %.4lf\n", NANOSECTOSEC(preprocessing_time_in_ns));
    // printf("Total time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    // printf("Memory cost (MB): %.4lf\n", BYTESTOMB(memory_cost_in_bytes));
    // printf("#Embeddings: %zu\n", embedding_count);
    // printf("Call Count: %zu\n", call_count);
    // printf("Per Call Count Time (nanoseconds): %.4lf\n", enumeration_time_in_ns / (call_count == 0 ? 1 : call_count));
    std::cout << "End." << std::endl;
    return 0;
}
