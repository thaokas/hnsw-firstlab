//
// Created by 35373 on 24-7-10.
//
#include "../../hnswlib/hnswlib.h"

int main(void)
{
    int _M_ = 32;
    int distribution_size = 400;
    size_t max_size_ = 2e8;
    size_t ramdom_seed = 100;
    double shaper = max_size_ * 1.0 / distribution_size;

    int Div_level_ = 0;
    while (shaper >= 1)
    {
        Div_level_ += 1;
        shaper /= _M_;
    }
    int Max_level_ = (log2(max_size_ * shaper) / log2(_M_)) + 1;


    std::cout << "shaper = " << shaper << " Max_level_ = " << Max_level_ <<  " Div_level_ = " << Div_level_ << std::endl;
    std::default_random_engine assit_lv_generator_;
    assit_lv_generator_.seed(ramdom_seed);
    int counter[20];
    memset(counter, 0, sizeof(counter));

    for (size_t i = 0; i < max_size_; ++i)
    {
        std::uniform_real_distribution<double> distribution(0.0, shaper);
        double r = -log2(distribution(assit_lv_generator_)) / log2(_M_);
        counter[int(r)] += 1;
    }


    for (auto item : counter)
        std::cout << "| " << item << " ";

    return 0;
}
