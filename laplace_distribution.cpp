#include "laplace_distribution.h"

LaplaceDistribution::LaplaceDistribution():
	n(1), mu(0), lambda(1){}

LaplaceDistribution::LaplaceDistribution(double _n, double _mu, double _lambda) :
	n(_n > 0 ? _n : throw 1), lambda(_lambda > 0 ? _lambda : throw 1), mu(_mu) {}

LaplaceDistribution::LaplaceDistribution(std::ifstream& file) {
	double n, mu, lambda;
	file.open("params.txt");
	if (!file.is_open()) {
		throw 0;
	}
	file >> n >> mu >> lambda;
	file.close();
	LaplaceDistribution(n, mu, lambda);
}

double LaplaceDistribution::standartization(const double x, const double lambda, const double mu) const {
	return (x - mu) / lambda;
}

double LaplaceDistribution::fact(const double n) const {
	if (n == 0 || n == 1) {
		return 1;
	}
	else {
		return n * fact(n - 1);
	}
}

double LaplaceDistribution::random_var() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

double LaplaceDistribution::density(const double x) const {
	double xstd = standartization(x, lambda, mu);
	double dist_sum = 0;
	for (int j = 0; j < n; ++j) {
		dist_sum += (fact(n - 1 + j) * pow(abs(xstd), (n - j - 1))) / (fact(n - j - 1) * fact(j) * pow(2, j));
	}
	return exp(-abs(xstd)) * dist_sum / (fact(n - 1) * pow(2, n) * lambda);
}

double LaplaceDistribution::expected_value() const {
	return mu;
}

double LaplaceDistribution::dispersion() const {
	return 2 * n * lambda * lambda;
}

double LaplaceDistribution::asymmetry() const {
	return 0;
}

double LaplaceDistribution::kurtosis() const {
	return 3 / n;
}

double LaplaceDistribution::rand_var() const {
	double mult1 = 1, mult2 = 1;
	for (int i = 0; i < n; ++i) {
		double r = random_var();
		if (r <= 0.5) {
			mult1 *= 2 * r;
		}
		if (r > 0.5) {
			mult2 *= 2 * (1 - r);
		}
	}
	return log(mult1 / mult2) * lambda + mu;
}

void LaplaceDistribution::set_form(const double n) {
	if (n <= 0) {
		throw 1;
	}
	this->n = n;
}

void LaplaceDistribution::set_shift(const double mu) {
	this->mu = mu;
}

void LaplaceDistribution::set_scale(const double lambda) {
	if (lambda <= 0) {
		throw 1;
	}
	this->lambda = lambda;
}

double LaplaceDistribution::get_form() const {
	return n;
}

double LaplaceDistribution::get_shift() const {
	return mu;
}

double LaplaceDistribution::get_scale() const {
	return lambda;
}

std::vector<double> LaplaceDistribution::generate_selection(const int n) const {
	std::vector<double> sample;
	for (int i = 0; i < n; ++i) {
		sample.push_back(rand_var());
	}
	sort(sample.begin(), sample.end());
	return sample;
}

std::vector<std::pair<double, double>> LaplaceDistribution::generate_graph_selection(const std::vector<double>& selection) const {
	std::vector<std::pair<double, double>> result;
	for (int i = 0; i < selection.size(); ++i) {
		result.push_back(std::make_pair(selection[i], density(selection[i])));
	}
	return result;
}

void LaplaceDistribution::save_in_file(std::ofstream& file) const {
	file.open("params.txt");
	file << n << std::endl << mu << std::endl << lambda;
	file.close();
}

void LaplaceDistribution::load_from_file(std::ifstream& file) {
	double n, mu, lambda;
	file.open("params.txt");
	if (!file.is_open()) {
		throw 0;
	}
	file >> n >> mu >> lambda;
	file.close();
	if (n <= 0 || lambda <= 0) {
		throw 1;
	}
	this->n = n;
	this->mu = mu;
	this->lambda = lambda;
}