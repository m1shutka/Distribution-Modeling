#pragma once
#include "distributions.h"

class LaplaceDistribution : public IDistribution, public IPresistend {
public:
	LaplaceDistribution();
	LaplaceDistribution(double n, double mu, double lambda);
	LaplaceDistribution(std::ifstream& file);

	double density(const double x) const override;
	double expected_value() const override;
	double dispersion() const override;
	double kurtosis() const override;
	double asymmetry() const override;
	double rand_var() const override;

	std::vector<double> generate_selection(const int n) const override;
	std::vector<std::pair<double, double>> generate_graph_selection(const std::vector<double>& selection) const override;

	void set_form(const double n);
	void set_shift(const double mu);
	void set_scale(const double lambda);

	double get_form() const;
	double get_shift() const;
	double get_scale() const;

	void load_from_file(std::ifstream& file) override;
	void save_in_file(std::ofstream& file) override;
private:
	double n;
	double mu;
	double lambda;

	double standartization(const double x, const double lambda, const double mu) const;
	double random_var() const;
	double fact(const double n) const;
};