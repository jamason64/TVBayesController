//dependencies:
#include <cmath>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <random> // for generation of data
#include <chrono>
#include <ctime>



using Eigen::MatrixXd;
using Eigen::VectorXd;

///retained stims
struct test_data{
    VectorXd time;
    MatrixXd inputs_samples;
    VectorXd output_samples;
};


///GPs
//multi_variable gaussian process
struct gaussian_process{
    MatrixXd domain;
    VectorXd mew;
    VectorXd sd;
};


//multi_variable gaussian process
struct gaussian_process_and_data{
    MatrixXd mew;
    MatrixXd sd;

};

//used for single variable model
struct gaussian_process_float{
    float mew;
    float sd;
};


///target stim
struct target_stim{
    VectorXd stim_variable;
    float expected_sig;
    float expected_mew;

};


///store change in tremor intensity
//float change_in_tremor_intensity;

///functions
float calc_kernel(VectorXd x1, VectorXd x2, float sigfsquared[3], float l[3]);
float get_result(VectorXd x, int time_iteration, std::mt19937 generate_noise, float noise_sd, VectorXd testvariables);
test_data initialise_data(int length_dataset, float x1_min, float x1_max, float x2_min, float x2_max, float x3_min, float x3_max, float noise_sd, VectorXd test_variables);
test_data uniform_time_increase(test_data dataset, float time_increase);
gaussian_process build_model(test_data dataset, int model_resolution_1, int model_resolution_2, int model_resolution_3,
                              float x1_min, float x1_max, float x2_min, float x2_max, float x3_min, float x3_max,
                              float model_decay_const);
target_stim aquisition_func(gaussian_process model, float aquisition_variable, int model_res_1, int model_res_2, int model_res_3, float random_number);
float evaluate_information(float expected_sig, float expected_mew, float true_change_trem_intensity, float decay_k, float K);
int replace_next_value(test_data dataset, gaussian_process model, target_stim aquisition_target, float change_in_tremor_intensity_);



