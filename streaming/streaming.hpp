#include "esa19.hpp"
#include <algorithm>
#include <string>
#include <unordered_set>
#include "../utils/dinic.hpp"
#include "../utils/ford_fulkerson.hpp"

// stats
double time_avg = 0;
double time_pivot_avg = 0;
double time_reps_avg = 0;
double time_matching_avg = 0;
double radius_eps_avg = 0;


// algorithm (Fair k-Center with Outliers)
class streaming {

    // dataset in instance
    std::vector<point> pset;

    // seed
    unsigned int seed = 0;

    // centers
    std::vector<unsigned int> center_set;

    // vector for dist-min
    std::vector<float> dist_min_array;

    // pivots
    std::vector<unsigned int> pivots;

    std::vector<std::vector<unsigned int>> reps;

    // group counter (in centers)
    std::unordered_map<std::string, unsigned int> group_counter;

    // group idx
    std::unordered_map<std::string, unsigned int> group_idx;
    std::unordered_map<unsigned int, std::string> group_idx_rev;

    // running time
    double running_time = 0;
    double running_time_pivot = 0;
    double running_time_reps = 0;
    double running_time_matching = 0;

    // radius
    float radius_guess = 0;
    float radius_eps = 0;

    
    // distance computation
    float compute_distance(const point *l, const point *r) {

        float distance = 0;
        for (unsigned int i = 0; i < dimensionality; ++i) distance += (l->pt[i] - r->pt[i]) * (l->pt[i] - r->pt[i]);
        return sqrt(distance);
    }

    // guessing OPT
    void guess_opt() {

        esa19 e(seed);
        radius_guess = e.fair_clustering(); // guess is done in consideration of outliers compared with no outlier case
    }

    // getPivots
    void get_pivots() {

        __start = std::chrono::system_clock::now();

        for (unsigned int i = 0; i < pset.size(); ++i) {

            bool f = 1;
            for (unsigned int j = 0; j < pivots.size(); ++j) {

                const float distance = compute_distance(&pset[i], &pset[pivots[j]]);
                if (distance <= 2 * radius_guess) {
                    f = 0;
                    break;
                }
            }
            if (f) {
                pivots.push_back(i);
                if (pivots.size() == k) break;
            }
        }

        __end = std::chrono::system_clock::now();
        running_time_pivot = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time_pivot /= 1000;
        running_time += running_time_pivot;
    }

    // getReps
    void get_reps() {

         __start = std::chrono::system_clock::now();

        reps.resize(pivots.size());
        std::vector<std::unordered_set<std::string>> reps_group(pivots.size());
        for (unsigned int i = 0; i < pivots.size(); ++i) {
            reps[i].push_back(pivots[i]);
            reps_group[i].insert(pset[pivots[i]].group);
        }

        for (unsigned int i = 0; i < pset.size(); ++i) {

            for (unsigned int j = 0; j < pivots.size(); ++j) {

                const float distance = compute_distance(&pset[i], &pset[pivots[j]]);
                if (distance <= radius_guess) {
                    if (reps_group[j].find(pset[i].group) == reps_group[j].end()) {
                        reps[j].push_back(i);
                        reps_group[j].insert(pset[i].group);
                    }
                }
            }
        }

        __end = std::chrono::system_clock::now();
        running_time_reps = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time_reps /= 1000;
        running_time += running_time_reps;
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

    // hittingSet
    void hitting_set() {

        __start = std::chrono::system_clock::now();

        get_group_idx();

        std::vector<std::vector<int>> graph(k + group_size + 2);
        for (unsigned int i = 0; i < k + group_size + 2; ++i) graph[i].resize(k + group_size + 2);
        for (unsigned int i = 1; i <= pivots.size(); ++i) graph[0][i] = 1;
        for (unsigned int i = k + 1; i < k + group_size + 1; ++i) graph[i][k + group_size + 1] = k_group_array[i - k - 1];
        for (unsigned int i = 0; i < reps.size(); ++i) {
            auto it = reps[i].begin();
            while (it != reps[i].end()) {
                graph[i+1][group_idx[pset[*it].group]] = 1;
                ++it;
            }
        }
        fordFulkerson( k + group_size + 2, graph, 0, k + group_size + 1);

        __end = std::chrono::system_clock::now();
        running_time_matching = std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();

        // make a flow graph
        Dinic<int> g(k + group_size + 2);
        for (unsigned int i = 0; i <= k; ++i) g.add_edge(0, i, 1);
        for (unsigned int i = k + 1; i < k + group_size + 1; ++i) g.add_edge(i, k + group_size + 1, k_group_array[i - k - 1]);
        for (unsigned int i = 0; i < reps.size(); ++i) {
            auto it = reps[i].begin();
            while (it != reps[i].end()) {
                g.add_edge(i+1, group_idx[pset[*it].group], 1);
                ++it;
            }
        }

        // get matching
        g.max_flow(0, k + group_size + 1);
        std::vector<std::pair<unsigned int, unsigned int>> match;
        g.match(match);

        __start = std::chrono::system_clock::now();

        for (unsigned int i = 0; i < match.size(); ++i) {
            const unsigned int idx = match[i].first - 1;
            
            // get match group
            std::string group_match = group_idx_rev[match[i].second];
            for (unsigned int j = 0; j < reps[i].size(); ++j) {
                if (pset[reps[i][j]].group == group_match) {
                    center_set.push_back(reps[i][j]);
                    
                    // update group counter
                    ++group_counter[group_match];
                    break;
                }
            }
        }

        __end = std::chrono::system_clock::now();
        running_time_matching += std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count();
        running_time_matching /= 1000;
        running_time += running_time_matching;
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

    // get remaining centers
    void get_remaining_centers() {

        __start = std::chrono::system_clock::now();

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

        __end = std::chrono::system_clock::now();
        running_time_pivot += (std::chrono::duration_cast<std::chrono::microseconds>(__end - __start).count()) / 1000;
        running_time += running_time_pivot;
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
    streaming() {}

    streaming(unsigned int s) {

        // init
        seed = s;
        pset = point_set;
        dist_min_array.resize(pset.size());
        const unsigned int size = pset.size();
        for (unsigned int i = 0; i < size; ++i) dist_min_array[i] = FLT_MAX;
    }

    // destructor
    ~streaming() {

        pset.shrink_to_fit();
        center_set.shrink_to_fit();
        dist_min_array.shrink_to_fit();
        pivots.shrink_to_fit();
        reps.shrink_to_fit();
    }

    // clustering
    void fair_clustering(const unsigned int seed) {

        // OPT guess
        guess_opt();

        // sort by random
        std::mt19937 engine(seed);
        std::shuffle(pset.begin(), pset.end(), engine);

        // get pivots
        get_pivots();

        // get reps
        get_reps();

        // hittingSet
        hitting_set();

        get_remaining_centers();
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
        
        f_name += "id(" + std::to_string(dataset_id) + ")_k(" + std::to_string(k) + ")_z(" + std::to_string(z) + ")_eps(" + std::to_string(eps) + ")_m(" + std::to_string(group_size) + ")_n(" + std::to_string(cardinality) + ")_streaming.csv";
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
            << ", run time for getPivot [msec]: " << running_time_pivot
            << ", run time for getReps [msec]: " << running_time_reps
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
        time_pivot_avg += running_time_pivot;
        time_reps_avg += running_time_reps;
        time_matching_avg += running_time_matching;
        radius_eps_avg += radius_eps;

        // averate result
        if (flag) file << "avg. run time [msec] " << time_avg / run_num
            << ", avg. run time (getPivot) [msec]: " << time_pivot_avg / run_num
            << ", avg. run time (getReps) [msec]: " << time_reps_avg / run_num
            << ", avg. run time (matching) [msec]: " << time_matching_avg / run_num
            << ", avg. radius-eps: " << radius_eps_avg / run_num
            << "\n";

        file.close();
    }
};