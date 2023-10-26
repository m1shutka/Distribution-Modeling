#include <iostream>
#include <fstream>
#include "laplace_distribution.h"
#include "empirical_distribution.h"
#include "mixture_distribution.cpp"
using namespace std;

int main() {
	setlocale(LC_ALL, "ru");
	ofstream file("save.txt");

	LaplaceDistribution ld1(1, 2, 3);
	LaplaceDistribution ld2(4, 5, 6);
	EmpiricalDistribution ed(ld1, 10, 1);
	MixtureDistribution<LaplaceDistribution, LaplaceDistribution> mx(ld1, ld2, 0.5);
	mx.component1().get_form();
	mx.density(0);
	
	MixtureDistribution<LaplaceDistribution, MixtureDistribution<LaplaceDistribution, LaplaceDistribution>> mx2(ld1, mx, 0.5);
	mx2.save_in_file(file);
	//file << typeid(mx2).name();
	file.close();

	ifstream file2("save.txt");

	mx2.load_from_file(file2);
	file2.close();
}