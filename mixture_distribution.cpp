#include "distributions.h"

template<class Distribution1, class Distribution2>
class MixtureDistribution : public IDistribution, public IPresistend {
public:
	MixtureDistribution(Distribution1& d1, Distribution2& d2, double p):
		d1(d1), d2(d2), p(p) {};

	double density(const double x) const override;
	double expected_value() const override;
	double dispersion() const override;
	double kurtosis() const override;
	double asymmetry() const override;
	double rand_var() const override;
	std::vector<double> generate_selection(const int n) const override;
	std::vector<std::pair<double, double>> generate_graph_selection(const std::vector<double>& selection) const override;

	void load_from_file(std::ifstream& file) override;
	void save_in_file(std::ofstream& file) override;

	double get_p() const;
	void set_p(const double p);

	Distribution1& component1() {return d1;}
	Distribution2& component2() {return d2;}
private:
	Distribution1 d1;
	Distribution2 d2;
	double p;

	double random_var() const;
};

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2>::get_p() const {
	return p;
}

template<class Dist1, class Dist2>
void MixtureDistribution<Dist1, Dist2>::set_p(const double p) {
	if (p < 0 or p > 1){
		throw 1;
	}
	this->p = p;
}

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2>::density(const double x) const {
	return (1 - p) * d1.density(x) + p * d2.density(x);
}

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2>::expected_value() const {
	return (1 - p) * d1.expected_value() + p * d2.expected_value();
}

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2>::dispersion() const {
	return (1 - p) * (pow(d1.expected_value(), 2) + d1.dispersion()) +
		p * (pow(d2.expected_value(), 2) + d2.dispersion()) -
		pow(expected_value(), 2);
}

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2>::asymmetry() const {
	return ((1 - p) * (pow((d1.expected_value() - expected_value()), 3) + 3 * (d1.expected_value() - expected_value()) * d1.dispersion() + pow(d1.dispersion(), 3 / 2) * d1.asymmetry()) +
		p * (pow((d2.expected_value() - expected_value()), 3) + 3 * (d2.expected_value() - expected_value()) * d2.dispersion() + pow(d2.dispersion(), 3 / 2) * d2.asymmetry())) /
		pow(dispersion(), 3 / 2);
}

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2> ::kurtosis() const {
	return ((1 - p) * (pow((d1.expected_value() - expected_value()), 4) + 6 * d1.dispersion() * pow((d1.expected_value() - expected_value()), 2) +
		4 * (d1.expected_value() - expected_value()) * pow(d1.dispersion(), 3 / 2) * d1.asymmetry() + pow(d1.dispersion(), 2) * d1.kurtosis()) +
		p * (pow((d2.expected_value() - expected_value()), 4) + 6 * d2.dispersion() * pow((d2.expected_value() - expected_value()), 2) +
			4 * (d2.expected_value() - expected_value()) * pow(d2.dispersion(), 3 / 2) * d2.asymmetry() + pow(d2.dispersion(), 2) * d2.kurtosis()) - 3) /
		pow(dispersion(), 2);
}

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2>::random_var() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

template<class Dist1, class Dist2>
double MixtureDistribution<Dist1, Dist2>::rand_var() const {
	double r = random_var();
	if (r > p) {
		return d1.rand_var();
	}
	else {
		return d2.rand_var();
	}
}

template<class Dist1, class Dist2>
std::vector<double> MixtureDistribution<Dist1, Dist2>::generate_selection(const int n) const {
	std::vector<double> sample;
	for (int i = 0; i < n; ++i) {
		sample.push_back(rand_var());
	}
	sort(sample.begin(), sample.end());
	return sample;
}

template<class Dist1, class Dist2>
std::vector<std::pair<double, double>> MixtureDistribution<Dist1, Dist2>::generate_graph_selection(const std::vector<double>& selection) const {
	std::vector<std::pair<double, double>> result;
	for (int i = 0; i < selection.size(); ++i) {
		result.push_back(std::make_pair(selection[i], density(selection[i])));
	}
	return result;
}

template<class Dist1, class Dist2>
void MixtureDistribution<Dist1, Dist2>::save_in_file(std::ofstream& file){
	file << p << std::endl;
	component1().save_in_file(file);
	component2().save_in_file(file);
	
}

template<class Dist1, class Dist2>
void MixtureDistribution<Dist1, Dist2>::load_from_file(std::ifstream& file) {
	file >> p;
	component1().load_from_file(file);
	component2().load_from_file(file);
}