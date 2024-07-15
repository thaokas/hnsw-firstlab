#include "../../hnswlib/hnswlib.h"
#include "hnswlib/assit.h"
#include <ctime>
#include <fstream>
#include <string>

void funtion(int _dim, int _max_elements, size_t _distribution_size) {
    int dim = _dim;               // Dimension of the elements
    int max_elements = _max_elements;   // Maximum number of elements, should be known beforehand
    int M = 32;                 // Tightly connected with internal dimensionality of the data
                                // strongly affects the memory consumption
    int ef_construction = 200;  // Controls index search speed/build speed tradeoff
    int query_size = max_elements / 20;
    hnswlib::divint div_strategy = 2;
    size_t distribution_size = _distribution_size;
    int k = distribution_size;

    // Initing index
    hnswlib::L2Space space(dim);
    hnswlib::HierarchicalNSW<float>* alg_hnsw = new hnswlib::HierarchicalNSW<float>(&space, max_elements, M, ef_construction, div_strategy, distribution_size);
    hnswlib::AlgorithmInterface<float>* alg_brute = new hnswlib::BruteforceSearch<float>(&space, max_elements);

    // Generate random data and query
    std::mt19937 rng;
    rng.seed(47);
    std::uniform_real_distribution<> distrib_real;
    float* data = new float[dim * max_elements];
    float* query = new float[dim * query_size];
    for (unsigned long long i = 0; i < dim * max_elements; ++i)
        data[i] = distrib_real(rng);
    for (unsigned long long i = 0; i < dim * query_size; ++i)
        query[i] = distrib_real(rng);

    // Add data to index
    clock_t start = clock();

    clock_t add_start = clock();
    for (hnswlib::labeltype i = 0; i < max_elements; i++)
    {
        alg_hnsw->addPoint(data + i * dim, i);
        alg_brute->addPoint(data + i * dim, i);
    }

    clock_t add_end = clock();
    std::cout << std::endl << "Time used for build database " << double(add_end - add_start) << std::endl;

    auto assist = new hnswlib::DatabaseStructure<float>();
    assist->print_db_structure(alg_hnsw);
    assist->print_db_size(alg_hnsw);

    // double correct = 0;
    clock_t search_start = clock();

    /*
     * level 1 400 点
     * level 1 1600个点 对 这些点聚类
     * ivf ？
     * 我的聚类中心
     */

    // Query the elements for themselves and measure recall
    std::vector<int> search_k_set;
    for (int i = 0; i < query_size; i++)
    {
        auto temp_result = alg_brute->searchKnnCloserFirst(query + i * dim, 1);
        auto query_res = alg_hnsw->searchInDivLayer(query + i * dim, k);
        hnswlib::labeltype gt_label = temp_result[0].second;
        auto gt_res = alg_hnsw->searchInDivLayer(data + gt_label * dim, 1);

        for (int i = 0; i < k; ++i)
        {
            if (query_res[i].second == gt_res[0].second)
            {
                search_k_set.push_back(i);
                // correct += 1;
            }
        }
    }

    clock_t search_end = clock();
    std::cout << "\nTime of search " << double(search_end - add_end) << std::endl;

    std::sort(search_k_set.begin(), search_k_set.end());

    std::cout << "Query size = " << search_k_set.size() << std::endl;
    std::cout << "Recall = 100% : " << search_k_set[search_k_set.size() - 1] << std::endl;
    std::cout << "Recall = 99.9% : " << search_k_set[(search_k_set.size() - 1) * 0.999] << std::endl;
    std::cout << "Recall = 99% : " << search_k_set[(search_k_set.size() - 1) * 0.99] << std::endl;

    // double recall = correct / query_size;
    // std::cout << "Recall: " << recall << "\n";


    // std::string out_file;
    // std::cin >> out_file;
    // std::ofstream result(out_file);
    // std::streambuf* cout_old = std::cout.rdbuf(result.rdbuf());
    //
    // assist->print_tag_info(alg_hnsw);
    //
    // std::cout.rdbuf(cout_old);

    delete[] data;
    delete[] query;
    delete alg_hnsw;
    delete alg_brute;
    // return 0;
}

int main(void)
{
    funtion(128, 200000, 400);
    funtion(512, 200000, 400);
    return 0;
}

