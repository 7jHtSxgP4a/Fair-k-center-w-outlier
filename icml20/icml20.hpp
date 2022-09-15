#include "../utils/data.hpp"
#include "../utils/dinic.hpp"
#include <functional>
#include <map>
#include <chrono>
#include <algorithm>
#include <functional>
#include <string>
#include <unordered_set>


// stats
double time_avg = 0;
double time_greedy_avg = 0;
double time_matching_avg = 0;
double radius_eps_avg = 0;
std::chrono::system_clock::time_point _s, _e;


// algorithm
class icml20 {

    // dataset in instance
    std::vector<point> pset;

    // seed
    unsigned int seed = 0;

    // centers
    std::vector<unsigned int> center_set;
    //std::vector<unsigned int> cluster;

    // vector for dist-min
    std::vector<float> dist_min_array;

    // group counter (in centers)
    std::unordered_map<std::string, unsigned int> group_counter;

    // vector for dist-to-center
    std::vector<float> dist_to_center;

    // group idx
    std::unordered_map<std::string, unsigned int> group_idx;
    std::unordered_map<unsigned int, std::string> group_idx_rev;


    // running time
    double running_time = 0;
    double running_time_greedy = 0;
    double running_time_matching = 0;

    // radius
    float radius_eps = 0;

   
    // distance computation
    float compute_distance(const point *l, const point *r) {

        float distance = 0;
        for (unsigned int i = 0; i < dimensionality; ++i) distance += (l->pt[i] - r->pt[i]) * (l->pt[i] - r->pt[i]);
        return sqrt(distance);
    }

    // gonzalez's algorithm
    bool gonzalez() {

        // get size
        const unsigned int size = pset.size();

        // prepare counter
        auto it = population.begin();
        while (it != population.end()) {
            group_counter.insert({it->first, 0});
            ++it;
        }

        // init dist-to-center
        dist_to_center.resize(k);
        dist_to_center[0] = FLT_MAX;

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

            float dist_max = 0;

            // update (1) dist-min to the intermediate result & (2) candidate
            for (unsigned int j = 0; j < size; ++j) {

                // (1) dist-min update
                const float distance = compute_distance(p, &pset[j]);
                if (dist_min_array[j] > distance) dist_min_array[j] = distance;

                // (2) furthest update
                if (dist_min_array[j] > dist_max) {

                    dist_max = dist_min_array[j];
                    idx = j;
                }
            }

            // determine next center
            center_set.push_back(idx);

            // update dist-to-center
            dist_to_center[i] = dist_max;

            // increment group counter
            ++group_counter[pset[idx].group];
        }

        // get last center
        point* p = &pset[idx];

        // update dist-min to the result
        for (unsigned int j = 0; j < size; ++j) {
            const float distance = compute_distance(p, &pset[j]);
            if (dist_min_array[j] > distance) dist_min_array[j] = distance;
        }

        /***********************/
        /*** violation check ***/
        /***********************/
        auto itr = k_group.begin();
        while (itr != k_group.end()) {
            if (group_counter[itr->first] > itr->second) return 1;
            ++itr;
        }

        return 0;
    }

    // get group index
    void get_group_idx() {

        // determine group idx
        unsigned int g_idx = k + 1;

        // transfomr k_group into array
        auto it = k_group.begin();
        while (it != k_group.end()) {

            k_group_array.push_back(it->second);
            group_idx[it->first] = g_idx;
            group_idx_rev[g_idx] = it->first;
            ++g_idx;
            ++it;
        }
    }

    // get largest h for shift
    unsigned int get_largest_h() {

        std::vector<std::pair<unsigned int, unsigned int>> match;

        // get size
        const unsigned int size = pset.size();

        unsigned int h_lo = 0;
        unsigned int h_hi = k;

        while (h_lo < h_hi) {

            unsigned int l = std::ceil((double)(h_lo + h_hi) / 2);

            Dinic<int> g(k + group_size + 2);
            for (unsigned int i = k + 1; i < k + group_size + 1; ++i) g.add_edge(i, k + group_size + 1, k_group_array[i - k - 1]);

            // update graph
            for (unsigned int j = 0; j < l; ++j) {

                g.add_edge(0, j + 1, 1);

                // compute dist-min to each group
                std::unordered_map<std::string, float> dist_min_group;
                auto it = k_group.begin();
                while (it != k_group.end()) {

                    dist_min_group[it->first] = FLT_MAX;
                    ++it;
                }

                for (unsigned int i = 0; i < size; ++i) {

                    // get group
                    const std::string group = pset[i].group;

                     // compute distance
                     const float distance = compute_distance(&pset[center_set[j]], &pset[i]);

                    // update dist-min to a given group
                    if (dist_min_group[group] > distance) dist_min_group[group] = distance;
                }

                it = k_group.begin();
                while (it != k_group.end()) {

                    if (dist_min_group[it->first] <= dist_to_center[l-1] / 2.0) g.add_edge(j + 1, group_idx[it->first], 1);
                    ++it;
                }
            }

            // max-flow
            int maxflow = g.max_flow(0, k + group_size + 1);
            if (maxflow == l) {
                //G_ = G;
                h_lo = l;
            }
            else {
                h_hi = l - 1;
            }
        }

        return h_hi;
    }

    // get remaining centers
    void get_remaining_centers(unsigned int idx) {

        if (center_set.size() < k) {

            // get size
            const unsigned int size = pset.size();

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

    // get match
    void get_match(const unsigned int h) {

        std::vector<std::pair<unsigned int, unsigned int>> match;

        // get size
        const unsigned int size = pset.size();

        // make a flow graph
        Dinic<int> g(k + group_size + 2);
        for (unsigned int i = k + 1; i < k + group_size + 1; ++i) g.add_edge(i, k + group_size + 1, k_group_array[i - k - 1]);
        for (unsigned int j = 0; j < h; ++j) {

            //G.addEdge(0, j + 1, 1);
            g.add_edge(0, j + 1, 1);

            // compute dist-min to each group
            std::unordered_map<std::string, float> dist_min_group;
            auto it = k_group.begin();
            while (it != k_group.end()) {

                dist_min_group[it->first] = FLT_MAX;
                ++it;
            }

            for (unsigned int i = 0; i < size; ++i) {

                // get group
                const std::string group = pset[i].group;

                // compute distance
                const float distance = compute_distance(&pset[center_set[j]], &pset[i]);

                // update dist-min to a given group
                if (dist_min_group[group] > distance) dist_min_group[group] = distance;
            }

            it = k_group.begin();
            while (it != k_group.end()) {

                if (dist_min_group[it->first] <= dist_to_center[h - 1] / 2.0) g.add_edge(j + 1, group_idx[it->first], 1);
                ++it;
            }
        }
        g.max_flow(0, k + group_size + 1);
        g.match(match);

        // clustering
        std::unordered_map<unsigned int, std::string> match_group;
        std::unordered_map<unsigned int, std::pair<float, unsigned int>> match_min;
        for (unsigned int i = 0; i < h; ++i) {
            match_group[match[i].first - 1] = group_idx_rev[match[i].second];
            match_min[match[i].first - 1] = {FLT_MAX, 0};
        }

        for (unsigned int i = 0; i < size; ++i) {

            // init dist_min_array
            dist_min_array[i] = FLT_MAX;

            double dist_min = FLT_MAX;
            unsigned int center_idx = 0;

            for (unsigned int j = 0; j < h; ++j) {

                const float distance = compute_distance(&pset[i], &pset[center_set[j]]);
                if (distance < dist_min) {
                    dist_min = distance;
                    center_idx = j;
                }
            }
            if (pset[i].group == match_group[center_idx]) {
                if (match_min[center_idx].first > dist_min) match_min[center_idx] = {dist_min, i};
            }
        }

        // clear center
        center_set.clear();

        // init group counter
        auto it = group_counter.begin();
        while (it != group_counter.end()) {
            it->second = 0;
            ++it;
        }

        auto itr = match_min.begin();
        while (itr != match_min.end()) {

            // get center
            center_set.push_back(itr->second.second);

            // increment group counter
            ++group_counter[pset[itr->second.second].group];
            ++itr;
        }

        // compute dist_min_array
        unsigned int idx = 0;
        float dist_min_max = 0;
        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned int j = 0; j < center_set.size(); ++j) {
                const float distance = compute_distance(&pset[i], &pset[center_set[j]]);
                if (distance < dist_min_array[i]) {
                    dist_min_array[i] = distance;
                    //cluster[i] = j;
                }
            }

            if (dist_min_max < dist_min_array[i]) {
                if (group_counter[pset[i].group] < k_group[pset[i].group]) {
                    dist_min_max = dist_min_array[i];
                    idx = i;
                }
            }
        }

        // get remaining centers
        get_remaining_centers(idx);
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
    icml20() {}

    icml20(unsigned int s) {

        // init
        seed = s;
        pset = point_set;
        dist_min_array.resize(pset.size());
        //cluster.resize(pset.size());

        const unsigned int size = pset.size();
        for (unsigned int i = 0; i < size; ++i) dist_min_array[i] = FLT_MAX;
    }

    // destructor
    ~icml20() {

        pset.shrink_to_fit();
        center_set.shrink_to_fit();
        dist_min_array.shrink_to_fit();
        dist_to_center.shrink_to_fit();
    }

    // clustering
    void fair_clustering() {

        _s = std::chrono::system_clock::now();

        // vanilla k-center clustering
        bool f = gonzalez();

        _e = std::chrono::system_clock::now();
        running_time_greedy = std::chrono::duration_cast<std::chrono::microseconds>(_e - _s).count();
        running_time_greedy /= 1000;
        running_time += running_time_greedy;

        _s = std::chrono::system_clock::now();

        // matching
        if (f) {

            // determine group idx
            get_group_idx();

            const unsigned int h = get_largest_h();
            get_match(h);
        }

        _e = std::chrono::system_clock::now();
        running_time_matching = std::chrono::duration_cast<std::chrono::microseconds>(_e - _s).count();
        running_time_matching /= 1000;
        running_time += running_time_matching;
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
        
        f_name += "id(" + std::to_string(dataset_id) + ")_k(" + std::to_string(k) + ")_z(" + std::to_string(z) + ")_eps(" + std::to_string(eps) + ")_m(" + std::to_string(group_size) + ")_n(" + std::to_string(cardinality) + ")_icml20.csv";
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

        // averate result
        if (flag) file << "avg. run time [msec] " << time_avg / run_num
            << ", avg. run time (greedy) [msec]: " << time_greedy_avg / run_num
            << ", avg. run time (matching) [msec]: " << time_matching_avg / run_num
            << ", avg. radius-eps: " << radius_eps_avg / run_num
            << "\n";
        file.close();
    }
};