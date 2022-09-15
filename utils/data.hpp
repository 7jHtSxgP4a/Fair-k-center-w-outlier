#include "file_input.hpp"
#include <cmath>
#include <cfloat>
#include <random>
#include <string>


// definition of a point
class point {

public:

    // identifier
    unsigned int identifier = 0;

    // vector
    std::vector<float> pt;

    // group
    std::string group = "";

    // flag (inlier or outlier)
    bool in = 1;

    // constructor
    point() {}

    point(const unsigned int id, std::vector<float> &p, std::string &g) {
        identifier = id;
        pt = p;
        group = g;
    }
};

// dataset
std::vector<point> point_set;
std::unordered_map<std::string, unsigned int> population;

// input dataset
void input_dataset() {

    // id variable
	unsigned int id = 0;

	// point coordinates variable
	point p;
    p.pt.resize(dimensionality);

    // min & max
    std::vector<float> min(dimensionality), max(dimensionality);
    for (unsigned int i = 0; i < dimensionality; ++i) {
        min[i] = FLT_MAX;
        max[i] = FLT_MIN;
    }

    // dataset input
	std::string f_name = "../_dataset/";
	if (dataset_id == 0) f_name += "adult_gender.csv";
    if (dataset_id == 1) f_name += "adult_race.csv";
    if (dataset_id == 2) f_name += "covertype.csv";
    if (dataset_id == 3) f_name += "diabete_gender.csv";
    if (dataset_id == 4) f_name += "diabete_race.csv";
    if (dataset_id == 5) f_name += "mirai.csv";
    if (dataset_id == 6) f_name += "kdd.csv";

    // file input
	std::ifstream ifs_file(f_name);
    std::string full_data;

	// error check
	if (ifs_file.fail()) {
		std::cout << " data file does not exist." << std::endl;
		std::exit(0);
	}

    std::mt19937 mt_(1);
    std::uniform_int_distribution<> rnd_m(0, syn_m - 1);

    // read data
	while (std::getline(ifs_file, full_data)) {

		std::string meta_info;
		std::istringstream stream(full_data);
		std::string type = "";

		for (unsigned int i = 0; i <= dimensionality; ++i) {

			std::getline(stream, meta_info, ',');
			std::istringstream stream_(meta_info);

            if (i < dimensionality) {
                const float val = (float)std::stold(meta_info);
                p.pt[i] = val;

                if (val < min[i]) min[i] = val;
                if (val > max[i]) max[i] = val;
            }
            else {
                if (dataset_id == 6 || dataset_id == 7) {
                    p.group = std::to_string(rnd_m(mt_));
                }
                else {
                    meta_info.pop_back();
                    p.group = meta_info;
                }
                
                if (population.find(p.group) == population.end()) {
                    population[p.group] = 1;
                }
                else {
                    population[p.group] += 1;
                }
            }
		}

		// update id
		p.identifier = id;

		// insert into dataset
		point_set.push_back(p);

		// increment identifier
		++id;
	}

    // shift max
    for (unsigned int i = 0; i < dimensionality; ++i) max[i] -= min[i];

    // get cardinality
    const unsigned int size = (unsigned int)point_set.size();

    // normalization
    const float max_coordinate = 100;
    for (unsigned int i = 0; i < size; ++i) {
        for (unsigned int j = 0; j < dimensionality; ++j) {
            point_set[i].pt[j] -= min[j];
            point_set[i].pt[j] /= max[j];
            point_set[i].pt[j] *= max_coordinate;
        }
    }

    // assign k_i + group size
    auto it = population.begin();
    while (it != population.end()) {
        ++group_size;
        unsigned int k_i = std::ceil(((float)it->second / size) * k);
        k_group[it->first] = k_i;
        ++it;
    }

    // random generator
    std::mt19937 mt(0);
	std::uniform_real_distribution<> rnd(-2 * max_coordinate, 3 * max_coordinate);
    std::uniform_int_distribution<> rnd_group(0, size - 1);

    // outlier injection
    p.in = 0;
    for (unsigned i = 0; i < z; ++i) {

        // determine coordinate
        while (1) {

            for (unsigned int j = 0; j < dimensionality; ++j) p.pt[j] = rnd(mt);
            bool f = 1;
            for (unsigned int j = 0; j < dimensionality; ++j) {
                if (p.pt[j] < 0 || p.pt[j] > max_coordinate) {
                    f = 0;
                    break;
                }
            }
            if (f == 0) break;
        }

        // determine group
        p.group = point_set[rnd_group(mt)].group;

        point_set.push_back(p);
    }


    // get cardinality
    cardinality = point_set.size() - z;

    std::cout << " ------------------\n";
	std::cout << " data id: " << dataset_id << "\n";
	std::cout << " dimensionality: " << dimensionality << "\n";
	std::cout << " cardinality: " << cardinality << "\n";
	std::cout << " k: " << k << "\n";
    it = population.begin();
    while (it != population.end()) {
        std::cout << " k for " << it->first << ": " << k_group[it->first] << " (" << it->first << " has " << it->second << " points)\n";
        ++it;
    }
	std::cout << " z: " << z << "\n";
    std::cout << " m: " << group_size << "\n";
    std::cout << " epsilon: " << eps << "\n";
	std::cout << " ------------------\n\n";
}