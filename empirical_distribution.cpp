#include "empirical_distribution.h"

EmpiricalDistribution::EmpiricalDistribution(const IDistribution& d, int n, int _k):
	size(n > 1 ? n : throw 1), k(_k > 2 ? _k : (int)log2(size) + 1) {
	selection = d.generate_selection(size);
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

EmpiricalDistribution::EmpiricalDistribution(const EmpiricalDistribution& ed):
	size(ed.size > 1 ? ed.size : throw 1), k(ed.k > 2 ? ed.k : ((int)log2(size) + 1)), selection(ed.selection), empirical_density(ed.empirical_density) {
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

EmpiricalDistribution::EmpiricalDistribution(std::ifstream& file) {
	file.open("selectionin.txt");
	if (!file.is_open()) {
		throw 0;
	}
	double x;
	while (file >> x) {
		selection.push_back(x);
	}
	size = selection.size();
	k = (int)log2(size) + 1;
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

EmpiricalDistribution::EmpiricalDistribution(const std::vector<double>& selection):
	size(selection.size()), k((int)log2(size) + 1), selection(selection) {
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

EmpiricalDistribution::~EmpiricalDistribution() {
	selection.clear();
	empirical_density.clear();
}

std::vector<double> EmpiricalDistribution::create_empirical_density(const double delta, std::vector<double> intervals) const {
	std::vector<double> result;
	int j = 0;
	for (int i = 0; i < intervals.size() - 1; ++i) {
		int count = 0;
		for (j; j < size; ++j) {
			if (i + 1 == intervals.size() - 1) {
				if (selection[j] <= intervals[i + 1]) {
					++count;
				}
				else {
					break;
				}
			}
			else {
				if (selection[j] < intervals[i + 1]) {
					++count;
				}
				else {
					break;
				}
			}
		}
		result.push_back(count / (size * delta));
	}
	return result;
}

double EmpiricalDistribution::random_var() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

double EmpiricalDistribution::qumulative_probability(const int i) const {
	double q = 0;
	for (int j = 0; j <= i; ++j) {
		q += empirical_density[j];
	}
	return q;
}

double EmpiricalDistribution::top_bound() const {
	double sum = 0;
	for (int i = 0; i < k; ++i) {
		sum += empirical_density[i];
	}
	return sum;
}

double EmpiricalDistribution::rand_var() const {
	double r;
	do r = (double)rand() / RAND_MAX * top_bound();
	while (r == 0 || r == top_bound());
	for (int i = 0; i < k - 1; ++i) {
		if (r > qumulative_probability(i) and r < qumulative_probability(i + 1)) {
			do r = (double)rand() / RAND_MAX * ((selection[0] + delta_calc() * (i + 1)) - (selection[0] + delta_calc() * i)) + (selection[0] + delta_calc() * (i + 1));
			while (r == (selection[0] + delta_calc() * i) || r == (selection[0] + delta_calc() * (i + 1)));
			break;
		}
	}
	return r;
}

EmpiricalDistribution& EmpiricalDistribution::operator = (const EmpiricalDistribution& ed) {
	if (this == &ed) {
		return *this;
	}
	selection = ed.selection;
	empirical_density = ed.empirical_density;
	size = ed.size;
	k = ed.k;
	return *this;
}

void EmpiricalDistribution::set_intrevals_number(const int k) {
	if (k >= 2) {
		this->k = k;
	}
	else {
		this->k = (int)log2(size) + 1;
	}
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

std::vector<double> EmpiricalDistribution::get_selection() const {
	return selection;
}

std::vector<double> EmpiricalDistribution::get_empirical_density() const {
	return empirical_density;
}

int EmpiricalDistribution::get_size() const {
	return size;
}

int EmpiricalDistribution::get_intrevals_number() const {
	return k;
}

double EmpiricalDistribution::delta_calc() const {
	return (selection[size - 1] - selection[0]) / (k);
}

std::vector<double> EmpiricalDistribution::create_intervals(const double delta) const {
	std::vector<double> intervals;
	double slider = selection[0];
	for (int i = 0; i < k + 1; ++i) {
		intervals.push_back(slider);
		slider += delta;
	}
	return intervals;
}

double EmpiricalDistribution::density(const double x) const {
	double slider = selection[0];
	if (x < slider) {
		return 0;
	}
	for (int i = 0; i < k; ++i) {
		slider += delta_calc();
		if (i == k - 1) {
			if (x <= slider) {
				return empirical_density[i];
			}
		}
		else {
			if (x < slider) {
				return empirical_density[i];
			}
		}
	}
	return 0;
}

double EmpiricalDistribution::expected_value() const {
	double sum = 0;
	for (int i = 0; i < size; ++i) {
		sum += selection[i];
	}
	return sum / size;
}

double EmpiricalDistribution::dispersion() const {
	double sum = 0;
	double exp_val = expected_value();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 2);
	}
	return sum / size;
}

double EmpiricalDistribution::asymmetry() const {
	double sum = 0;
	double exp_val = expected_value();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 3);
	}
	return sum / (size * pow(dispersion(), 3 / 2));
}

double EmpiricalDistribution::kurtosis() const {
	double sum = 0;
	double exp_val = expected_value();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 4);
	}
	return (sum / (size * pow(dispersion(), 2))) - 3;
}

std::vector<double>  EmpiricalDistribution::generate_selection(const int n) const {
	std::vector<double> result;
	for (int i = 0; i < n; ++i) {
		result.push_back(rand_var());
	}
	sort(result.begin(), result.end());
	return result;
}

std::vector<std::pair<double, double>> EmpiricalDistribution::generate_graph_selection(const std::vector<double>& selection) const {
	std::vector <std::pair<double, double >> result;
	for (int j = 0; j < size; ++j) {
		result.push_back(std::make_pair(selection[j], density(selection[j])));
	}
	return result;
}

void EmpiricalDistribution::save_in_file(std::ofstream& file){
	file.open("epmric_distribution.txt");
	for (int i = 0; i < size; ++i) {
		file << selection[i] << " " << std::endl;
	}
	file.close();
}

void EmpiricalDistribution::load_from_file(std::ifstream& file) {
	file.open("emparams.txt");
	if (!file.is_open()) {
		throw 0;
	}
}