#include "esa19.hpp"
#include <cmath>
#include <string>
#include <unordered_set>
#include "../utils/dinic.hpp"


// stats
double time_avg = 0;
double time_greedy_avg = 0;
double time_matching_avg = 0;
double radius_eps_avg = 0;
double radius_sqr = 0;

// algorithm (Fair k-Center with Outliers)
class fkco {

    // dataset in instance
    std::vector<point> pset;
    std::vector<unsigned int> sset;

    // seed
    unsigned int seed = 0;

    // centers
    std::vector<unsigned int> center_set;

    // vector for dist-min
    std::vector<float> dist_min_array;

    // group counter (in centers)
    std::unordered_map<std::string, unsigned int> group_counter;

    // group idx
    std::unordered_map<std::string, unsigned int> group_idx;
    std::unordered_map<unsigned int, std::string> group_idx_rev;

    std::vector<std::unordered_map<std::string, std::pair<float, unsigned int>>> group_in_center;

    // running time
    double running_time = 0;
    double running_time_greedy = 0;
    double running_time_matching = 0;

    // radius
    float radius_guess = 0;
    float radius_eps = 0;

    // sample size
    unsigned int sample_size = 0;

    
    // distance computation
    float compute_distance(const point *l, const point *r) {

        float distance = 0;
        for (unsigned int i = 0; i < dimensionality; ++i) distance += (l->pt[i] - r->pt[i]) * (l->pt[i] - r->pt[i]);
        return sqrt(distance);
    }

    // get clustering result
    unsigned int get_clustering_result() {

        // get size
        const unsigned int size = pset.size();

        unsigned int idx = 0;
        float dist_min_max = 0;
        for (unsigned int i = 0; i < size; ++i) {

            // init dist_min
            dist_min_array[i] = FLT_MAX;

            for (unsigned int j = 0; j < center_set.size(); ++j) {
                const float distance = compute_distance(&pset[i], &pset[center_set[j]]);
                if (distance < dist_min_array[i]) {
                    dist_min_array[i] = distance;
                }    
            }

            if (dist_min_max < dist_min_array[i]) {
                if (group_counter[pset[i].group] < k_group[pset[i].group]) {
                    dist_min_max = dist_min_array[i];
                    idx = i;
                }
            }
        }

        return idx;
    }

    // greedy iteration
    bool greedy() {

        bool flag = 1;

        __start = std::chrono::system_clock::now();

        // prepare counter
        auto it = population.begin();
        while (it != population.end()) {
            group_counter.insert({it->first, 0});
            ++it;
        }

        /******************/
        /*** 1st sample ***/
        /******************/
        std::mt19937 mt(seed);
        std::uniform_int_distribution<> rnd_int(0, sample_size - 1);
        unsigned int idx = rnd_int(mt);
        center_set.push_back(sset[idx]);

        // increment group counter
        ++group_counter[pset[sset[idx]].group];

        /**********************/
        /*** i-th iteration ***/
        /**********************/
        for (unsigned int i = 1; i <= k; ++i) {

            std::vector<unsigned int> candidate;

            // get last center
            point* p = &pset[sset[idx]];

            // update (1) dist-min to the intermediate result & (2) candidate
            for (unsigned int j = 0; j < sample_size; ++j) {

                // (1) dist-min update
                const float distance = compute_distance(p, &pset[sset[j]]);
                if (dist_min_array[j] > distance) dist_min_array[j] = distance;

                // (2) candidate update
                if (dist_min_array[j] > 2.0 * radius_guess && i < k) candidate.push_back(j);
            }

            const unsigned int candidate_size = candidate.size();
            if (candidate_size == 0) break;

            // sampling center
            std::uniform_int_distribution<> rnd_int_sample(0, candidate_size - 1);
            idx = candidate[rnd_int_sample(mt)];

            // update center
            center_set.push_back(sset[idx]);

            // increment group counter
            ++group_counter[pset[sset[idx]].group];
        }

        if (center_set.size() < k) {
            flag = 0;
        }
        else {

            /***********************/
            /*** violation check ***/
            /***********************/
            auto itr = k_group.begin();
            while (itr != k_group.end()) {
                if (group_counter[itr->first] > itr->second) {
                    flag = 0;
                    break;
                }
                ++itr;
            }
        }

        group_in_center.resize(center_set.size());
        const unsigned int size_ = pset.size();
        for (unsigned int i = 0; i < size_; ++i) {
            std::pair<float, unsigned int> to_assign = {FLT_MAX, 0};

            for (unsigned int j = 0; j < center_set.size(); ++j) {

                const float distance = compute_distance(&pset[i], &pset[center_set[j]]);
                if (to_assign.first > distance) to_assign = {distance, j};
                if (to_assign.first <= radius_guess) break;
            }

            if (to_assign.first <= radius_guess) {
                if (group_in_center[to_assign.second].find(pset[i].group) == group_in_center[to_assign.second].end()) {
                    group_in_center[to_assign.second].insert({pset[i].group, {to_assign.first, i}});
                }
            }
        }

        __end = std::chrono::system_clock::now();
        running_time_greedy = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time_greedy /= 1000;
        running_time += running_time_greedy;

        return flag;
    }

    // get group index
    void get_group_idx() {

        // determine group idx
        unsigned int g_idx = k + 1;

        // transform k_group into array
        auto it = k_group.begin();
        while (it != k_group.end()) {

            k_group_array.push_back(it->second);
            group_idx[it->first] = g_idx;
            group_idx_rev[g_idx] = it->first;
            ++g_idx;
            ++it;
        }
    }

    // get remaining centers
    void get_remaining_centers() {

        if (center_set.size() < k) {

            // get size
            const unsigned int size = pset.size();

            unsigned int idx = get_clustering_result();

            center_set.push_back(idx);

            // increment group counter
            ++group_counter[pset[idx].group];

            while (center_set.size() < k) {

                // get last center
                point* p = &pset[idx];

                // init dist-max
                float dist_min_max = 0;

                // update (1) dist-min to the intermediate result & (2) candidate
                for (unsigned int j = 0; j < size; ++j) {

                    // (1) dist-min update
                    const float distance = compute_distance(p, &pset[j]);
                    if (dist_min_array[j] > distance) dist_min_array[j] = distance;

                    // (2) candidate update
                    if (dist_min_max < dist_min_array[j]) {
                        if (group_counter[pset[j].group] < k_group[pset[j].group]) {
                            dist_min_max = dist_min_array[j];
                            idx = j;
                        }
                    }
                }

                center_set.push_back(idx);

                // increment group counter
                ++group_counter[pset[idx].group];
            }
        }
    }

    // get matching
    void matching() {

        __start = std::chrono::system_clock::now();

        /****************/
        /*** matching ***/
        /****************/
        get_group_idx();

        // make a flow graph
        Dinic<int> g(k + group_size + 2);
        for (unsigned int i = 0; i <= k; ++i) g.add_edge(0, i, 1);
        for (unsigned int i = k + 1; i < k + group_size + 1; ++i) g.add_edge(i, k + group_size + 1, k_group_array[i - k - 1]);
        for (unsigned int i = 0; i < group_in_center.size(); ++i) {
            auto it = group_in_center[i].begin();
            while (it != group_in_center[i].end()) {
                g.add_edge(i+1, group_idx[it->first], 1);
                ++it;
            }
        }

        __end = std::chrono::system_clock::now();
        double temp = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time_greedy += temp / 1000;
        __start = std::chrono::system_clock::now();

        // get matching
        g.max_flow(0, k + group_size + 1);
        std::vector<std::pair<unsigned int, unsigned int>> match;
        g.match(match);

        for (unsigned int i = 0; i < match.size(); ++i) {
            const unsigned int idx = match[i].first - 1;
            
            // get match group
            std::string group_match = group_idx_rev[match[i].second];
            if (pset[center_set[idx]].group != group_match) {

                // update group counter
                --group_counter[pset[center_set[idx]].group];

                // shift
                center_set[idx] = group_in_center[idx][group_match].second;

                // update group counter
                ++group_counter[group_match];
            }
        }

        __end = std::chrono::system_clock::now();
        running_time_matching = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time_matching /= 1000;
        running_time += running_time_matching;
        __start = std::chrono::system_clock::now();

        // supplement remaining centers
        get_remaining_centers();

        __end = std::chrono::system_clock::now();
        temp = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time += temp / 1000;
    }

    // get result stats
    void get_stats() {

        /**********************/
        /*** get max radius ***/
        /**********************/
        get_clustering_result();
        std::sort(dist_min_array.begin(), dist_min_array.end(), std::greater<float>());
        unsigned int idx = (unsigned int)((1.0 + eps) * z);
        radius_eps = dist_min_array[idx];
    }


public:

    // constructor
    fkco() {}

    fkco(unsigned int s, const float r_guess) {

        // init
        seed = s;
        pset = point_set;
        dist_min_array.resize(pset.size());
        const unsigned int size = pset.size();
        for (unsigned int i = 0; i < size; ++i) dist_min_array[i] = FLT_MAX;

        radius_guess = r_guess / 2.0;

        // sampling
        sample_size = (std::log(size) * size * k) / ((1 + eps) * (1 + eps) * z);

        std::mt19937 mt(seed * 10);
        std::uniform_int_distribution<> rnd_int(0, size - 1);
        std::unordered_set<unsigned int> test;
        while (test.size() < sample_size) {
            test.insert(rnd_int(mt));
        }
        auto it = test.begin();
        while (it != test.end()) {
            sset.push_back(*it);
            ++it;
        }
    }

    // destructor
    ~fkco() {

        pset.shrink_to_fit();
        sset.shrink_to_fit();
        center_set.shrink_to_fit();
        dist_min_array.shrink_to_fit();
    }

    // clustering
    void fair_clustering() {

        // greedy iteration
        bool flag = greedy();

        // matching for fairness
        if (flag == 0) matching();
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
        if (dataset_id == 6) f_name += "6_kdd/";
        
        f_name += "id(" + std::to_string(dataset_id) + ")_k(" + std::to_string(k) + ")_z(" + std::to_string(z) + ")_eps(" + std::to_string(eps) + ")_m(" + std::to_string(group_size) + ")_n(" + std::to_string(cardinality) + ")_ours_sample.csv";
        std::ofstream file;
        file.open(f_name.c_str(), std::ios::out | std::ios::app);

        if (file.fail()) {
            std::cerr << " cannot open the output file." << std::endl;
            file.clear();
            return;
        }

        // get stats
        get_stats();

        file << "run time [msec]: " << running_time
            << ", run time for greedy [msec]: " << running_time_greedy
            << ", run time for matching [msec]: " << running_time_matching
            << ", radius-eps: " << radius_eps
            << ", #samples: " << sample_size
            << ", k: " << center_set.size();
        auto it = group_counter.begin();
        while (it != group_counter.end()) {
            file << ", " << it->first << ": " << it->second;
            ++it;
        }
        file << "\n";

        time_avg += running_time;
        time_greedy_avg += running_time_greedy;
        time_matching_avg += running_time_matching;
        radius_eps_avg += radius_eps;
        radius_sqr += radius_eps * radius_eps;

        // averate result
        if (flag) {
            double temp = radius_sqr - (radius_eps_avg * radius_eps_avg / run_num);
            file << "avg. run time [msec] " << time_avg / run_num
            << ", avg. run time (greedy) [msec]: " << time_greedy_avg / run_num
            << ", avg. run time (matching) [msec]: " << time_matching_avg / run_num
            << ", avg. radius-eps: " << radius_eps_avg / run_num
            << ", std dev.: " << temp / (run_num - 1)
            << "\n";
        }
        file.close();
    }
};