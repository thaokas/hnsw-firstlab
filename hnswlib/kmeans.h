#pragma once

namespace hnswlib {
    template<typename dist_t>
    class KMeamsSearch : public AlgorithmInterface<dist_t>
    {
    public:
        char *data_;
        size_t maxelements_;
        size_t cur_element_count;
        size_t size_per_element_;
        int k_centroids;

        size_t data_size_;
        DISTFUNC <dist_t> fstdistfunc_;
        void *dist_func_param_;
        std::mutex index_lock;

        std::unordered_map<labeltype, size_t> dict_external_to_internal;

        KMeamsSearch(SpaceInterface<dist_t> *s)
            : data_(nullptr), maxelements_(0), cur_element_count(0), size_per_element_(0),
            data_size_(0), dist_func_param_(nullptr) {}

        KMeamsSearch(SpaceInterface <dist_t> *s, size_t maxElements, int centroids_number) {
            maxelements_ = maxElements;
            k_centroids = centroids_number;
            data_size_ = s->get_data_size();
            fstdistfunc_ = s->get_dist_func();
            dist_func_param_ = s->get_dist_func_param();
            size_per_element_ = data_size_ + sizeof(labeltype);
            data_ = (char *) malloc(maxElements * size_per_element_);
            if (data_ == nullptr)
                throw std::runtime_error("Not enough memory: KMeansSearch failed to allocate data");
            cur_element_count = 0;
        }

        ~KMeamsSearch()
        {
            free(data_);
        }

        void addPoint(const void *datapoint, labeltype label, bool replace_deleted = false)
        {
            int idx;
            {
                std::unique_lock<std::mutex> lock(index_lock);

                auto search = dict_external_to_internal.find(label);
                if (search != dict_external_to_internal.end()) {
                    idx = search->second;
                } else {
                    if (cur_element_count >= maxelements_) {
                        throw std::runtime_error("The number of elements exceeds the specified limit\n");
                    }
                    idx = cur_element_count;
                    dict_external_to_internal[label] = idx;
                    cur_element_count++;
                }
            }
            memcpy(data_ + size_per_element_ * idx + data_size_, &label, sizeof(labeltype));
            memcpy(data_ + size_per_element_ * idx, datapoint, data_size_);
        }

    };
}
