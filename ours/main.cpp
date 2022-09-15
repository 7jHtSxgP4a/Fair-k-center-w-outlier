//#include "ours.hpp"
#include "ours_sampling.hpp"


int main() {

    // file input
    input_parameter();

    // data input
    input_dataset();

    // display current time
    get_current_time();

    // random generator
    std::mt19937 mt(0);
	std::uniform_int_distribution<> rnd(0, 1.0);

    /******************/
    /*** guessing r ***/
    /******************/
    float r_guess = FLT_MAX;
    const unsigned int pre_num = 20;
    for (unsigned int i = 0; i < pre_num; ++i) {
        esa19 e(i);
        const float r = e.fair_clustering();
        if (r_guess > r) r_guess = r;
    }

    /*************************/
    /*** run our algorithm ***/
    /*************************/
    for (unsigned int i = 0; i < run_num; ++i) {

        // flag for last result
        bool flag = 0;
        if (i == run_num - 1) flag = 1;

        // make an instance
        fkco f(i, r_guess);

        // run
        f.fair_clustering();

        // result output
        f.output_file(flag);
    }

    return 0;
}