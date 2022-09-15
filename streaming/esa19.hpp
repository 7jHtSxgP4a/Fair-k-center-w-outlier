#include "../utils/data.hpp"
#include <functional>
#include <map>
#include <chrono>
#include <algorithm>
#include <functional>


// stats
std::chrono::system_clock::time_point __start, __end;


// algorithm for OPT guess
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
    float radius = 0;

    // distance computation
    float compute_distance(const point *l, const point *r) {

        float distance = 0;
        for (unsigned int i = 0; i < dimensionality; ++i) distance += (l->pt[i] - r->pt[i]) * (l->pt[i] - r->pt[i]);
        return sqrt(distance);
    }

    // candidate update
    void update_candidate(std::multimap<float, unsigned int, std::greater<float>> &candidate, float distance, unsigned idx) {

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


public:
    // constructor
    esa19() {}

    esa19(unsigned int s) {

        // init
        seed = s;
        pset = point_set;
        for (unsigned int i = 0; i < z; ++i) pset.pop_back();
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

    // get run time
    double get_run_time() { return running_time; }

    // clustering
    float fair_clustering() {

        __start = std::chrono::system_clock::now();

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
            const unsigned int candidate_size = candidate.size();
            std::uniform_int_distribution<> rnd_int_sample(0, candidate_size - 1);
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

        __end = std::chrono::system_clock::now();
        running_time = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time /= 1000;


        /*****************************/
        /*** get clsutering result ***/
        /*****************************/

        // get last center
        point* p = &pset[idx];

        // update dist-min to the result
        for (unsigned int j = 0; j < size; ++j) {
            const float distance = compute_distance(p, &pset[j]);
            if (dist_min_array[j] > distance) dist_min_array[j] = distance;
        }

        std::sort(dist_min_array.begin(), dist_min_array.end(), std::greater<float>());
        idx = (unsigned int)((1.0 + eps) * z);
        radius = dist_min_array[idx];

        return radius;
    }
};