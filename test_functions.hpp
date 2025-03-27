#include <cmath>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <random> // for generation of data
#include <chrono>
#include <ctime>


using Eigen::MatrixXd;
using Eigen::VectorXd;


struct ranges{
    float phase_m[2];
    float phase_a[2];
    float phase_b[2];

    float freq_a[2];
    float freq_b[2];
    float freq_s[2];

    float pw_a[2];
    float pw_b[2];
    float pw_s[2];
};


VectorXd get_test_var(VectorXd previous, ranges dist_range, int testoption, std::mt19937 gen_rd, float stepsize);
