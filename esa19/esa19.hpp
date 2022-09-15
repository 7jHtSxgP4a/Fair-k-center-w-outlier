#include "../utils/data.hpp"
#include <functional>
#include <map>
#include <chrono>
#include <algorithm>
#include <functional>
#include <string>


// stats
double time_avg = 0;
double radius_eps_avg = 0;
std::chrono::system_clock::time_point start, end;


// algorithm
class esa19 {

    // dataset in instance
    std::vector<point> pset;

    // seed
    unsigned int seed = 0;

    // centers
    std::vector<unsigned int> center_set;

    // vector for dist-min
    std::vector<float> dist_min_array;

    // group counter (in centers)
    std::unordered_map<std::string, unsigned int> group_counter;


    // running time
    double running_time = 0;

    // radius
    float radius_eps = 0;

    
    // distance computation
    float compute_distance(const point *l, const point *r) {

        float distance = 0;
        for (unsigned int i = 0; i < dimensionality; ++i) distance += (l->pt[i] - r->pt[i]) * (l->pt[i] - r->pt[i]);
        return sqrt(distance);
    }

    // candidate update
    void update_candidate(std::multimap<float, unsigned int, std::greater<float>> &candidate, const float distance, const unsigned idx) {

        const std::string group = pset[idx].group;
        const unsigned int cnt = group_counter[group];

        if (k_group[group] > cnt) {

            if (candidate.size() < (1.0 + eps) * z) {
                candidate.insert({distance,idx});
            }
            else {
                auto itr = candidate.end();
                --itr;
                if (itr->first < distance) {
                    candidate.insert({distance,idx});
                    itr = candidate.end();
                    candidate.erase(--itr);
                }
            }
        }
    }

    // get result stats
    void get_stats() {

        /**********************/
        /*** get max radius ***/
        /**********************/
        std::sort(dist_min_array.begin(), dist_min_array.end(), std::greater<float>());
        unsigned int idx = (unsigned int)((1.0 + eps) * z);
        radius_eps = dist_min_array[idx];
    }

public:
    // constructor
    esa19() {}

    esa19(unsigned int s) {

        // init
        seed = s;
        pset = point_set;
        dist_min_array.resize(pset.size());

        const unsigned int size = pset.size();
        for (unsigned int i = 0; i < size; ++i) dist_min_array[i] = FLT_MAX;
    }

    // destructor
    ~esa19() {

        pset.shrink_to_fit();
        center_set.shrink_to_fit();
        dist_min_array.shrink_to_fit();
    }

    // clustering
    void fair_clustering() {

        start = std::chrono::system_clock::now();

        // get size
        const unsigned int size = pset.size();

        // prepare counter
        auto it = population.begin();
        while (it != population.end()) {
            group_counter.insert({it->first, 0});
            ++it;
        }

        // random generator
        std::mt19937 mt(seed);
        std::uniform_int_distribution<> rnd_int(0, size - 1);

        /******************/
        /*** 1st sample ***/
        /******************/
        unsigned int idx = rnd_int(mt);
        center_set.push_back(idx);

        // increment group counter
        ++group_counter[pset[idx].group];

        /**********************/
        /*** i-th iteration ***/
        /**********************/
        for (unsigned int i = 1; i < k; ++i) {

            std::multimap<float, unsigned int, std::greater<float>> candidate;

            // get last center
            point* p = &pset[idx];

            // update (1) dist-min to the intermediate result & (2) candidate
            for (unsigned int j = 0; j < size; ++j) {

                // (1) dist-min update
                const float distance = compute_distance(p, &pset[j]);
                if (dist_min_array[j] > distance) dist_min_array[j] = distance;

                // (2) candidate update
                update_candidate(candidate, dist_min_array[j], j);
            }

            // sampling
            std::uniform_int_distribution<> rnd_int_sample(0, z + group_size);  // aim at sampling a furthest point for each group
            int rand_num = rnd_int_sample(mt);

            auto itr = candidate.begin();
            while (rand_num > 0) {
                ++itr;
                --rand_num;
            }

            idx = itr->second;
            center_set.push_back(idx);

            // increment group counter
            ++group_counter[pset[idx].group];
        }

        end = std::chrono::system_clock::now();
        running_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        running_time /= 1000;

        /*****************************/
        /*** get clustering result ***/
        /*****************************/

        // get last center
        point* p = &pset[idx];

        // update dist-min to the result
        for (unsigned int j = 0; j < size; ++j) {
            const float distance = compute_distance(p, &pset[j]);
            if (dist_min_array[j] > distance) dist_min_array[j] = distance;
        }
    }

    // output result
    void output_file(bool flag) {

        std::string f_name = "result/";
        if (dataset_id == 0) f_name += "0_adult-gender/";
        if (dataset_id == 1) f_name += "1_adult-race/";
        if (dataset_id == 2) f_name += "2_covertype/";
        if (dataset_id == 3) f_name += "3_diabetes-gender/";
        if (dataset_id == 4) f_name += "4_diabetes-race/";
        if (dataset_id == 5) f_name += "5_mirai/";
        if (dataset_id == 6) f_name += "5_kdd/";
        
        f_name += "id(" + std::to_string(dataset_id) + ")_k(" + std::to_string(k) + ")_z(" + std::to_string(z) + ")_eps(" + std::to_string(eps) + ")_m(" + std::to_string(group_size) + ")_n(" + std::to_string(cardinality) + ")_esa19.csv";
        std::ofstream file;
        file.open(f_name.c_str(), std::ios::out | std::ios::app);

        if (file.fail()) {
            std::cerr << " cannot open the output file." << std::endl;
            file.clear();
            return;
        }

        // get stats
        get_stats();

        file << "run time [msec]: " << running_time << ", radius_eps: " << radius_eps;
        auto it = group_counter.begin();
        while (it != group_counter.end()) {
            file << ", " << it->first << ": " << it->second;
            ++it;
        }
        file << "\n";

        time_avg += running_time;
        radius_eps_avg += radius_eps;

        // averate result
        if (flag) file << "avg. run time [msec]: " << time_avg / run_num
            << ", avg. radius-eps: " << radius_eps_avg / run_num
            << "\n";

        file.close();
    }
};