#pragma once
#include "distributions.h"

class EmpiricalDistribution : public IDistribution, public IPresistend {
public:
	EmpiricalDistribution(const IDistribution& d, int n, int k);
	EmpiricalDistribution(const EmpiricalDistribution& d);
	EmpiricalDistribution(std::ifstream& file);
	EmpiricalDistribution(const std::vector<double>& selection);
	~EmpiricalDistribution();

	double density(const double x) const override;
	double expected_value() const override;
	double dispersion() const override;
	double kurtosis() const override;
	double asymmetry() const override;
	double rand_var() const override;

	std::vector<double> generate_selection(const int n) const override;
	std::vector<std::pair<double, double>> generate_graph_selection(const std::vector<double>& selection) const override;

	EmpiricalDistribution& operator=(const EmpiricalDistribution& ed);

	void set_intrevals_number(int k);

	std::vector<double> get_selection() const;
	std::vector<double> get_empirical_density() const;
	int get_size() const;
	int get_intrevals_number() const;

	void load_from_file(std::ifstream& file) override;
	void save_in_file(std::ofstream& file) const override;
private:
	std::vector<double> selection;
	std::vector<double> empirical_density;
	int size;
	int k;

	double delta_calc() const;
	std::vector<double> create_intervals(const double delta) const;
	double random_var() const;
	std::vector<double> create_empirical_density(const double delta, std::vector<double> intervals) const;
	double qumulative_probability(const int i) const;
	double top_bound() const;
};