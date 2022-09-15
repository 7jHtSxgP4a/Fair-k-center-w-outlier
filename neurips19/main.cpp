#include "neurips19.hpp"


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

    // iterate run_num
    for (unsigned int i = 0; i < run_num; ++i) {

        // flag for last result
        bool f = 0;
        if (i == run_num - 1) f = 1;

        // make an instance
        neurips19 n(i);

        // run
        n.fair_clustering();

        // result output
        n.output_file(f);
    }


    return 0;
}