// This is a test file for testing the interface
//  >>> virtual std::vector<std::pair<dist_t, labeltype>>
//  >>>    searchKnnCloserFirst(const void* query_data, size_t k) const;
// of class AlgorithmInterface

#include "../../hnswlib/hnswlib.h"

#include <assert.h>

#include <vector>
#include <iostream>

namespace {

using idx_t = hnswlib::labeltype;

void test() {
    int d = 16;
    idx_t n = 10000;
    idx_t nq = 100;
    size_t k = 5;
    size_t max_elements = 2 * n;
    hnswlib::divint lv_div_strategy = 0;
    size_t device_number = 30;

    std::vector<float> data(n * d);
    std::vector<float> query(nq * d);

    std::mt19937 rng;
    rng.seed(47);
    std::uniform_real_distribution<> distrib;

    for (idx_t i = 0; i < n * d; ++i) {
        data[i] = distrib(rng);
    }
    for (idx_t i = 0; i < nq * d; ++i) {
        query[i] = distrib(rng);
    }

    hnswlib::L2Space space(d);
    hnswlib::AlgorithmInterface<float>* alg_brute  = new hnswlib::BruteforceSearch<float>(&space, max_elements);
    hnswlib::AlgorithmInterface<float>* alg_hnsw = new hnswlib::HierarchicalNSW<float>(&space, max_elements, 16, 200, lv_div_strategy, device_number);

    for (size_t i = 0; i < n; ++i) {
        alg_brute->addPoint(data.data() + d * i, i);
        alg_hnsw->addPoint(data.data() + d * i, i);
    }

    for (size_t j = 0; j < nq; ++j) {
        const void* p = query.data() + j * d;
        auto res_h = alg_hnsw->searchKnn(p, k);
        auto res_b = alg_brute->searchKnn(p, k);
        assert(res_b.size() == res_h.size());
        size_t t = res_h.size();
        for (int i = 0; i < t; ++i)
        {
            assert(res_b.top() == res_h.top());
            std::cout << "[ " << res_h.top().first << " | " << res_b.top().first << " ]" << std::endl;
            res_h.pop();
            res_b.pop();
        }
    }

    delete alg_brute;
    delete alg_hnsw;
}

}  // namespace

int main() {
    std::cout << "Testing ..." << std::endl;
    test();
    std::cout << "Test ok" << std::endl;

    return 0;
}
