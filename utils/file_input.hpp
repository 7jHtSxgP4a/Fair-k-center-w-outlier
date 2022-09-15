#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
#include <unordered_map>


// dataset identifier
unsigned int dataset_id = 0;

// data dimensionality
unsigned int dimensionality = 2;

// cluster size
unsigned int k = 2;
std::vector<unsigned int> k_group_array;	// set of k_i
std::unordered_map<std::string, unsigned int> k_group;	// set of k_i

// outlier number
float z = 0;

// group size
unsigned int group_size = 0;
unsigned int syn_m = 0;

// cardinality for synthetic data
unsigned int cardinality = 0;

// outlier extension
float eps = 0;

// run number
const unsigned int run_num = 100;

// get current time
void get_current_time() {

	time_t t = time(NULL);
	printf(" %s\n\n", ctime(&t));
}

// parameter input
void input_parameter() {

	std::ifstream ifs_dataset_id("parameter/dataset_id.txt");
	std::ifstream ifs_k("parameter/k.txt");
	std::ifstream ifs_z("parameter/z.txt");
	std::ifstream ifs_eps("parameter/epsilon.txt");
	std::ifstream ifs_m("parameter/m.txt");
	std::ifstream ifs_n("parameter/n.txt");

	if (ifs_dataset_id.fail()) {
		std::cout << " dataset_id.txt does not exist." << std::endl;
		std::exit(0);
	}
	if (ifs_k.fail()) {
		std::cout << " k.txt does not exist." << std::endl;
		std::exit(0);
	}
	if (ifs_z.fail()) {
		std::cout << " z.txt does not exist." << std::endl;
		std::exit(0);
	}
	if (ifs_eps.fail()) {
		std::cout << " epsilon.txt does not exist." << std::endl;
		std::exit(0);
	}
	if (ifs_m.fail()) {
		std::cout << " m.txt does not exist." << std::endl;
		std::exit(0);
	}
	if (ifs_n.fail()) {
		std::cout << " n.txt does not exist." << std::endl;
		std::exit(0);
	}

	while (!ifs_dataset_id.eof()) { ifs_dataset_id >> dataset_id; }
	while (!ifs_k.eof()) { ifs_k >> k; }
	while (!ifs_z.eof()) { ifs_z >> z; }
	while (!ifs_eps.eof()) { ifs_eps >> eps; }
	while (!ifs_m.eof()) { ifs_m >> syn_m; }
	while (!ifs_n.eof()) { ifs_n >> cardinality; }

	// set dimensionality
	if (dataset_id == 0 || dataset_id == 1) dimensionality = 6;
	if (dataset_id == 2) dimensionality = 10;
	if (dataset_id == 3 || dataset_id == 4) dimensionality = 7;
	if (dataset_id == 5) dimensionality = 115;
	if (dataset_id == 6) dimensionality = 16;
}