#include "streaming.hpp"


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


    /*********************/
    /*** run streaming ***/
    /*********************/
    for (unsigned int i = 0; i < run_num; ++i) {

        // flag for last result
        bool flag = 0;
        if (i == run_num - 1) flag = 1;

        // make an instance
        streaming s(i);

        // run
        s.fair_clustering(i);

        // result output
        s.output_file(flag);
    }

    return 0;
}