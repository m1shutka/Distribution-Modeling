#define CATCH_CONFIG_RUNNER
#include <iostream>
#include <fstream>
#include "laplace_distribution.h"
#include "empirical_distribution.h"
#include "mixture_distribution.cpp"
#include "catch.hpp"
using namespace std;

int run_unit_tests(int argc, char** argv) {
	int result = Catch::Session().run(argc, argv);
	return result;
}

int main(int argc, char** argv) {
	setlocale(LC_ALL, "ru");
	srand((unsigned)time(0));
	ofstream file1("sdt_distr.txt");
	ofstream file2("emp_distr.txt");
	
	try {
		/// ����� ���
		LaplaceDistribution ld1(1, -6, 2);
		LaplaceDistribution ld2(1, -2, 1);
		LaplaceDistribution ld3(1, 2, 1);
		LaplaceDistribution ld4(1, 6, 2);
		MixtureDistribution<LaplaceDistribution, LaplaceDistribution> md1(ld1, ld2, 0.5);
		MixtureDistribution<LaplaceDistribution, LaplaceDistribution> md2(ld3, ld4, 0.5);
		MixtureDistribution<MixtureDistribution<LaplaceDistribution, LaplaceDistribution>, MixtureDistribution<LaplaceDistribution, LaplaceDistribution>> md(md2, md1, 0.5);
		auto selection = md.generate_selection(1000);
		auto graph1 = md.generate_graph_selection(selection);
		EmpiricalDistribution ed(selection);
		auto graph2 = ed.generate_graph_selection(selection);
		EmpiricalDistribution ed2(ed, 1000, 1);
		auto graph3 = ed2.generate_graph_selection(selection);

		cout << "���������: n1 = " << md.component1().component1().get_form() << ", mu1 = " << md.component1().component1().get_shift() << ", lambda1 = " << md.component1().component1().get_scale()<<
			"\nn2 = " << md.component1().component2().get_form() << ", mu2 = " << md.component1().component2().get_shift() << ", lambda2 = " << md.component1().component2().get_scale() <<
			"\nn3 = " << md.component2().component1().get_form() << ", mu3 = " << md.component2().component1().get_shift() << ", lambda3 = " << md.component2().component1().get_scale() <<
			"\nn4 = " << md.component2().component2().get_form() << ", mu4 = " << md.component2().component2().get_shift() << ", lambda4 = " << md.component2().component2().get_scale() << ", p = " << md.get_p() << endl << endl;;
		cout << "������������� ��������������:\n���. ��������:" << md.expected_value() << "\n���������: " << md.dispersion() << "\n����. ����������: " << md.asymmetry() << "\n����. ��������:" << md.kurtosis() << endl << endl;
		cout << "������������ ��������������1:\n���. ��������:" << ed.expected_value() << "\n���������: " << ed.dispersion() << "\n����. ����������: " << ed.asymmetry() << "\n����. ��������:" << ed.kurtosis() << endl << endl;
		cout << "������������ ��������������2:\n���. ��������:" << ed2.expected_value() << "\n���������: " << ed2.dispersion() << "\n����. ����������: " << ed2.asymmetry() << "\n����. ��������:" << ed2.kurtosis() << endl;



		for (int i = 0; i < 1000; ++i) {
			file1 << graph2[i].first << "\t" << graph2[i].second << endl;
			file2 << graph3[i].first << "\t" << graph3[i].second << endl;
		}
		///
	}
	catch (const int error) {
		if (error == 0) {
			cout << "��� ����� ��� ������� ���������� � ������ params.txt!!!" << endl;
		}
		else if (error == 1) {
			cout << "����������� ������� ��������� �������������!!!" << endl;
		}
	}
	file1.close();
	file2.close();
	run_unit_tests(argc, argv);
}