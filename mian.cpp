#include <iostream>
#include "laplace_distribution.h"
#include "empirical_distribution.h"
#include "mixture_distribution.cpp"
using namespace std;

int main() {
	setlocale(LC_ALL, "ru");

	LaplaceDistribution ld1(1, 0, 1);
	LaplaceDistribution ld2(1, 0, 1);
	EmpiricalDistribution ed(ld1, 10, 1);
	MixtureDistribution<LaplaceDistribution, LaplaceDistribution> mx(ld1, ld2, 0.5);
	mx.component1().get_form();
}