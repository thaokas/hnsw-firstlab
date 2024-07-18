#pragma once

namespace hnswlib {
    typedef size_t centroidtype;

    template<typename dist_t>
    class KMeamsSearch : public AlgorithmInterface<dist_t>
    {
    public:
        char *data_;
        char* centroids_data_;
        size_t maxelements_;
        size_t cur_element_count;
        size_t size_per_element_;
        size_t size_per_centroid_;
        int k_centroids;

        size_t data_size_;
        DISTFUNC <dist_t> fstdistfunc_;
        void *dist_func_param_;
        std::mutex index_lock;

        std::unordered_map<labeltype, size_t> dict_external_to_internal;
        std::unordered_map<labeltype, centroidtype> district_lool_up_;

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
            size_per_centroid_ = data_size_;
            data_ = (char *) malloc(maxElements * size_per_element_);
            if (data_ == nullptr)
                throw std::runtime_error("Not enough memory: KMeansSearch failed to allocate data");
            centroids_data_ = (char *) malloc(k_centroids * size_per_centroid_);
            if (centroids_data_ == nullptr)
                throw std::runtime_error("Not enough memory: KMeansSearch failed to allocate centroid");
            cur_element_count = 0;
            std::cout << "Init end\n";
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

                // search ： 根据外部标签确定内部索引
                auto search = dict_external_to_internal.find(label);
                if (search != dict_external_to_internal.end()) {
                    // 更新向量
                    idx = search->second;
                } else {
                    if (cur_element_count >= maxelements_)
                        throw std::runtime_error("The number of elements exceeds the specified limit\n");

                    idx = cur_element_count;
                    dict_external_to_internal[label] = idx;
                    cur_element_count++;
                }
            }
            memcpy(data_ + size_per_element_ * idx + data_size_, &label, sizeof(labeltype));
            memcpy(data_ + size_per_element_ * idx, datapoint, data_size_);
        }

        char * getCentroidDataByIndex(centroidtype index) const
        {
            return centroids_data_ + index * size_per_centroid_;
        }

        char * getDataByIndex(int idx) const
        {
            return data_ + idx * size_per_element_;
        }

        labeltype getExternalLabelByIndex(int idx)
        {
            return *(labeltype *)(data_ + size_per_element_ * idx + data_size_);
        }

        void generateKMeans()
        {
            std::default_random_engine rng;
            rng.seed(32);
            std::uniform_real_distribution<float> distrib_real;
            float *src_cent_ptr = (float *)centroids_data_;
            for (centroidtype i = 0; i < k_centroids * data_size_ / sizeof(float); ++i)
            {
                src_cent_ptr[i] = distrib_real(rng);
            }

            std::unordered_map<centroidtype, std::vector<size_t>> means_n_times;

            bool changed = true;
            int rec_times = 0;
            while (changed)
            {
                rec_times += 1;
                for (size_t i = 0; i < cur_element_count; ++i)
                {
                    char * currObj = data_ + size_per_element_ * i;
                    centroidtype minIndex = 0;
                    dist_t d = fstdistfunc_(currObj, getCentroidDataByIndex(minIndex), dist_func_param_);
                    for (centroidtype j = 1; j < k_centroids; ++j)
                    {
                        dist_t temp_d = fstdistfunc_(currObj, getCentroidDataByIndex(j), dist_func_param_);
                        if (d > temp_d)
                        {
                            d = temp_d;
                            minIndex = j;
                        }
                    }
                    means_n_times[minIndex].push_back(i);
                }

                char * tmp_centroid_data = (char *) malloc(k_centroids * size_per_centroid_);
                memset(tmp_centroid_data, 0, sizeof(tmp_centroid_data));

                for (centroidtype i = 0; i < k_centroids; ++i)
                {
                    for (int j = 0; j < means_n_times[i].size(); ++i)
                    {
                        float * tptr =(float*)(tmp_centroid_data + size_per_centroid_ * i);
                        float * tdata = (float*)getDataByIndex(means_n_times.find(i)->second[j]);
                        for (size_t k = 0; k < *(size_t*)dist_func_param_; ++k)
                        {
                            tptr[k] = (tptr[k] * j + tdata[k]) / (j + 1);
                        }
                    }
                }

                // for (centroidtype i = 0; i < k_centroids; ++i)
                // {
                //     std::cout << "[ centroid id = " << i << " | " << means_n_times[i].size() << " ]\n";
                // }

                if (memcmp(tmp_centroid_data, centroids_data_, k_centroids * data_size_) == 0)
                {
                    changed = false;
                }
                else
                {
                    delete[] centroids_data_;
                    centroids_data_ = tmp_centroid_data;
                    means_n_times.clear();
                }
            }

            for (centroidtype i = 0; i < k_centroids; ++i)
            {
                for (auto item : means_n_times[i])
                {
                    district_lool_up_[getExternalLabelByIndex(item)] = i;
                }

                std::cout << "[ centroid id = " << i << " | " << means_n_times[i].size() << " ]\n";
            }
            std::cout << "------" << rec_times << " times-------\n";
        }

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchKnn(const void*, size_t, BaseFilterFunctor* isIdAllowed = nullptr) const
        {
            std::cout << "..." << std::endl;
            std::priority_queue<std::pair<dist_t, labeltype>> res;
            return res;
        }


        void saveIndex(const std::string &location)
        {
            std::cout << "..." << std::endl;
        }
    };
}
