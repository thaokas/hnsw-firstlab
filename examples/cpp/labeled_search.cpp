//
// Created by 35373 on 24-7-12.
//

#include "../../hnswlib/hnswlib.h"
#include "hnswlib/assit.h"
#include <ctime>
#include <fstream>
#include <string>
#include <cassert>

int main() {
    // 将标准输出重定位到文件输出


    int dim = 16;               // Dimension of the elements
    int max_elements = 200000;   // Maximum number of elements, should be known beforehand
    int M = 32;                 // Tightly connected with internal dimensionality of the data
                                // strongly affects the memory consumption
    int ef_construction = 40;  // Controls index search speed/build speed tradeoff
    int query_size = max_elements / 20;

    // Initing index
    hnswlib::L2Space space(dim);
    hnswlib::HierarchicalNSW<float>* alg_hnsw = new hnswlib::HierarchicalNSW<float>(&space, max_elements, M, ef_construction, 2, 5);
    hnswlib::AlgorithmInterface<float>* alg_brute  = new hnswlib::BruteforceSearch<float>(&space, max_elements);

    // Generate random data
    std::mt19937 rng;
    rng.seed(47);
    std::uniform_real_distribution<> distrib_real;
    float* data = new float[dim * max_elements];
    for (unsigned long long i = 0; i < dim * max_elements; i++) {
        data[i] = distrib_real(rng);
    }

    // Add data to index
    clock_t start = clock();
    for (int i = 0; i < max_elements; i++) {
        alg_hnsw->addPoint(data + i * dim, i);
        alg_brute->addPoint(data + i * dim, i);
    }
    clock_t add_end = clock();

    auto assist = new hnswlib::DatabaseStructure<float>();
    assist->print_db_structure(alg_hnsw);
    assist->print_db_size(alg_hnsw);

    std::cout << std::endl << "Time used for build database " << double(add_end - start) << std::endl;

    alg_hnsw->NearestLabeler();
    clock_t label_end = clock();

    std::cout << std::endl << "Time used for labeling " << double(label_end - add_end) << std::endl;

    std::vector<int> devices_recorder;
    size_t sum_for_avg = 0;

    // generate query data
    float *query = new float[query_size * dim];
    for (unsigned long long i = 0; i < query_size * dim; ++i)
    {
        query[i] = distrib_real(rng);
    }

    // Query the elements for themselves and measure recall
    double correct = 0;
    for (int i = 0; i < query_size; i++) {
        auto true_result = alg_brute->searchKnn(query + i * dim, 1);
        auto query_res = alg_hnsw->searchNLayerST(true_result.top().second, query + i * dim, alg_hnsw->ef_, alg_hnsw->level_div_);
        // auto gt_result = alg_hnsw->searchNLayerST(,)
    }
    double recall = correct / query_size;
    std::cout << "Recall: " << recall << "\n";

    clock_t search_end = clock();

    std::cout << "\nTime of search " << double(search_end - label_end) << std::endl;

    std::cout << "\nSearch " << sum_for_avg * 1.0 / devices_recorder.size() << " devices in average" << std::endl << std::endl;
    // for (auto item : devices_recorder)
    // {
    //     std::cout << item << " | ";
    // }

    delete[] data;
    delete alg_hnsw;
    return 0;
}
