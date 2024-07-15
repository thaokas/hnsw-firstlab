#pragma once
#include <iomanip>

namespace hnswlib
{
    template <typename T>
    class DatabaseStructure
    {
     public:

        void print_db_structure(HierarchicalNSW<T>* alg_hnsw)
        {
            std::cout << "db cur_element_count: " << alg_hnsw->cur_element_count << std::endl;
            std::vector<int>  each_level_count(alg_hnsw->maxlevel_ + 1, 0);
            std::cout << "Below are info about levels" << std::endl;

            for (int i = 0; i < alg_hnsw->cur_element_count; ++i)
                each_level_count[alg_hnsw->element_levels_[i]] += 1;
            std::cout << std::endl;

            for (int i = 0; i <= alg_hnsw->maxlevel_; ++i)
                std::cout << "There are " << each_level_count[i] << " vectors start from level " << i << std::endl;
        }

        void print_db_size(HierarchicalNSW<T>* alg_hnsw)
        {
            std::cout << "| M_ = " << alg_hnsw->M_ << " | maxlevel = " << alg_hnsw->maxlevel_ << " |" << std::endl;
            unsigned long long sum_bytes = alg_hnsw->maxlevel_ * alg_hnsw->size_links_level0_;
            for (auto item : alg_hnsw->element_levels_)
            {
                sum_bytes += item * alg_hnsw->size_links_per_element_;
            }
            std::cout << "sum size of db = " << sum_bytes << " bytes" << std::endl;
        }

        void print_tag_info(HierarchicalNSW<T>* alg_hnsw)
        {
            std::unordered_map<labeltype, std::unordered_set<labeltype>> childs_record_;
            std::cout << std::endl << std::setw(30) << std::setfill('-') << "Each's Domainer" << std::setw(20) << std::setfill('-') << "-" << std::endl;
            for (auto item : alg_hnsw->domainer_lookup_)
            {
                if (item.second.size() > 1)
                    std::cout << item.second.size() << " | ";
                for (auto ele : item.second)
                {
                    childs_record_[ele].emplace(item.first);
                }
            }
            std::cout << std::endl << std::setw(30) << std::setfill('-') << "Domainer's Children" << std::setw(20) << std::setfill('-') << "-" << std::endl;
            std::cout << std::endl;
            size_t sum_extra = 0;
            for (auto item : childs_record_)
            {
                std::cout << item.second.size() << " | ";
                sum_extra += item.second.size();
            }
            std::cout << "\n" << "Sum size = " << sum_extra << std::endl;
        }
    };
}


