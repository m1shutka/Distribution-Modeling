#include "catch.hpp"
#include "laplace_distribution.h"
#include "mixture_distribution.cpp"
#include "empirical_distribution.h"

bool equal(const double& x, const double& y) {
    if (abs(x - y) <= 0.1) {
        return true;
    }
    else return false;
}

TEST_CASE("basic methods") {
    LaplaceDistribution distr;
    CHECK(distr.get_form() == 1);
    CHECK(distr.get_shift() == 0);
    CHECK(distr.get_scale() == 1);
}

TEST_CASE("standart distribution") {
    LaplaceDistribution distr;
    CHECK(distr.density(0) == 0.5);
    CHECK(distr.expected_value() == 0);
    CHECK(distr.dispersion() == 2);
    CHECK(distr.asymmetry() == 0);
    CHECK(distr.kurtosis() == 3);
}

TEST_CASE("shift scale transformation") {
    LaplaceDistribution distr;
    distr.set_scale(2);
    distr.set_shift(2);
    CHECK(equal(distr.density(0), 0.091) == true);
    CHECK(distr.expected_value() == 2);
    CHECK(distr.dispersion() == 8);
    CHECK(distr.asymmetry() == 0);
    CHECK(distr.kurtosis() == 3);
}

TEST_CASE("mixture distribution") {
    LaplaceDistribution distr1;
    LaplaceDistribution distr2;
    distr1.set_scale(2);
    distr2.set_scale(2);
    distr1.set_shift(2);
    distr2.set_shift(2);
    MixtureDistribution<LaplaceDistribution, LaplaceDistribution> mdistr(distr1, distr2, 0.5);
    
    CHECK(equal(mdistr.density(0), 0.091) == true);
    CHECK(mdistr.expected_value() == 2);
    CHECK(mdistr.dispersion() == 8);
    CHECK(mdistr.asymmetry() == 0);
    CHECK(equal(mdistr.kurtosis(), 2.953) == true);
}

TEST_CASE("mixture distribution expected") {
    LaplaceDistribution distr1;
    LaplaceDistribution distr2;
    distr1.set_scale(2);
    distr2.set_scale(2);
    distr1.set_shift(1);
    distr2.set_shift(2);
    MixtureDistribution<LaplaceDistribution, LaplaceDistribution> mdistr(distr1, distr2, 0.5);
    CHECK(mdistr.expected_value() == 1.5);
}

TEST_CASE("mixture distribution dispersion") {
    LaplaceDistribution distr1;
    LaplaceDistribution distr2;
    distr1.set_scale(1);
    distr2.set_scale(3);
    MixtureDistribution<LaplaceDistribution, LaplaceDistribution> mdistr(distr1, distr2, 0.5);
    CHECK(mdistr.dispersion() == 10);
}

TEST_CASE("late binding mechanism") {
    LaplaceDistribution distr1;
    LaplaceDistribution distr2;
    distr1.set_scale(1);
    distr2.set_scale(3);
    
    MixtureDistribution<LaplaceDistribution, LaplaceDistribution> mdistr(distr1, distr2, 0.5);
    MixtureDistribution<LaplaceDistribution, MixtureDistribution<LaplaceDistribution, LaplaceDistribution>> mdistr2(distr1, mdistr, 0.5);
    CHECK(mdistr2.component1().get_form()== 1);
    CHECK(mdistr2.component1().get_shift() == 0);
    CHECK(mdistr2.component1().get_scale() == 1);
    CHECK(mdistr2.component2().component1().get_form() == 1);
    CHECK(mdistr2.component2().component1().get_shift() == 0);
    CHECK(mdistr2.component2().component1().get_scale() == 1);
    CHECK(mdistr2.component2().component2().get_form() == 1);
    CHECK(mdistr2.component2().component2().get_shift() == 0);
    CHECK(mdistr2.component2().component2().get_scale() == 3);
}

TEST_CASE("empirical distribution") {
    LaplaceDistribution distr1;
    EmpiricalDistribution ed(distr1, 200, 1);
    CHECK(ed.get_size() == 200);
    CHECK(ed.get_intrevals_number() == 8);
    CHECK(equal(ed.expected_value(), 0.0) == true);
}