#pragma once
#include <vector>
#include <fstream>
#include <algorithm>

class IDistribution {
public:
	double virtual density(const double x) const = 0;
	double virtual expected_value() const = 0;
	double virtual dispersion() const = 0;
	double virtual kurtosis() const = 0;
	double virtual asymmetry() const = 0;
	double virtual rand_var() const = 0;
	std::vector<double> virtual generate_selection(const int n) const = 0;
	std::vector<std::pair<double, double>> virtual generate_graph_selection(const std::vector<double>& selection) const = 0;
};

class IPresistend {
	void virtual load_from_file(std::ifstream& file) = 0;
	void virtual save_in_file(std::ofstream& file) = 0;
};