#include "gplibrary.hpp"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::last;
using Eigen::all;

//This code is a C++ implementation of TV-GP-UCB bayesian optimisation adapted to adapt 
//  its exploration bias to reflect the difference between the expected information gain
//  and the real information again (exploration is increased if true uncertainty is more 
//  than the real uncertainty).
// For detailed information consult:
//  Applying Bayesian Optimisation as a Controller for Peripheral Stimulation,
//  James Mason,
//  University of Oxford, Honour School of Engineering Science, Part C Project,
//  2023
//Contact: james.w.mason@gmail.com

//usage (see toplevel for an implementation used in testing):

// ///initialisation of remembered samples feeding model (equivalent of in-clinic set-up)
// retained_dataset = initialise_data(remember_n_stims,
//                                     phase_rad_min, phase_rad_max,
//                                     freq_hz_min, freq_hz_max,
//                                         pulse_duration_min_us, pulse_duration_max_us,
//                                         standard_deviation_noise, test_var);
// for (int run_n = 0; run_n < number_of_optimisation_iterations; run_n++){
//     model = build_model(retained_dataset, model_resolution_phase, model_resolution_freq, model_resolution_pulse,
//                             phase_rad_min, phase_rad_max,
//                             freq_hz_min, freq_hz_max,
//                             pulse_duration_min_us, pulse_duration_max_us,
//                             model_time_decay);
//     next_stimulation = aquisition_func(model, aquisition_variable, model_resolution_phase, model_resolution_freq, model_resolution_pulse, distnv(gen));
//     change_in_tremor_intensity = get_result(next_stimulation.stim_variable, run_n, gen, standard_deviation_noise, test_var);
//     aquisition_variable = evaluate_information(next_stimulation.expected_sig, next_stimulation.expected_mew, change_in_tremor_intensity, decay_constant, aquisition_proportional_const);
//     replacement_index = replace_next_value(retained_dataset, model,
//                                             next_stimulation, change_in_tremor_intensity); 
//     retained_dataset = uniform_time_increase(retained_dataset, 1);
//     retained_dataset.inputs_samples(replacement_index,0) = next_stimulation.stim_variable[0];
//     retained_dataset.inputs_samples(replacement_index,1) = next_stimulation.stim_variable[1];
//     retained_dataset.inputs_samples(replacement_index,2) = next_stimulation.stim_variable[2];
//     retained_dataset.output_samples(replacement_index) = change_in_tremor_intensity;
//     retained_dataset.time(replacement_index) = 0;
// }



///calculate kernel
float calc_kernel(VectorXd x1, VectorXd x2, float sigfsquared[3], float l[3]){
    // kernel function between two variables
    // This calculates the difference between the two vectors
    // This can be done in any way but a good guide is: https://www.cs.toronto.edu/~duvenaud/cookbook/
    // Current implementation is to evaluate an elementwise distance and then combine for kernel output
    // input is the two comparison vectors (EIGEN library) and two float vectors which guide the length
    //       scale and output variance.
    // output is a scalar float

    float a;
    float b;
    float v[3];

    //phase:
    float x1temp = x1(0);
    float x2temp = x2(0);
    a = x1temp - x2temp;
    b = 0.5*abs(a);
    a = sin(b)/l[0];
    b = -2*pow(a,2);
    a = exp(b);
    v[0] = sigfsquared[0] *a;

    //voltage
    x1temp = x1(1);
    x2temp = x2(1);
    a = x1temp - x2temp;
    b = abs(a);
    a = (b)/(l[1]);
    b = -pow((a/2), 2);
    a = exp(b);
    v[1] = sigfsquared[1] * a;

    //pulse duration
    x1temp = x1(2);
    x2temp = x2(2);
    a = x1temp - x2temp;
    b = abs(a);
    a = (b)/(l[1]);
    b = -pow((a/2),2);
    a = exp(b);
    v[2] = sigfsquared[2] * a;


    float k = v[0]*(v[1]*v[2]);
    //cout<<"\n"<<v[0]<<"__"<<v[1]<<"__"<<v[2];
    return k;

    ///kernel functions:

    /// RBF kernel
    //VectorXd v = x1 - x2;
    //float a = v.squaredNorm();
    //float b = pow(a, 2);
    //a = -(b)/(2*l);
    //b = exp(a);
    //float k = sigfsquared * b;
    //return k;


    /// periodic kernel // see kernel cookbook duvenaud
    //VectorXd v = x1 - x2;
    //float a = v.squaredNorm();
    //float b = a/2;
    //a = sin(b)/l;
    //b = -2*pow(a,2);
    //a = exp(b);
    //float k = sigfsquared *a;
    //return k;

    ///time based kernel (RBF atm)
    //VectorXd v = x1 - x2;
    //float a = v.squaredNorm();
    //float b = pow(a, 2);
    //a = -(b)/(2*l);
    //b = exp(a);
    //float k = sigfsquared * b;
    //return k;
}


///perform stimulation, return change in tremor intensity - This would be changed out for an actual result acquisition 
float get_result(VectorXd x, int time_iteration, std::mt19937 generate_noise, float noise_sd, VectorXd testvariables){
    //get_result is a place holder for data acqusition
    //atm this outputs a scalar float for a one dimensional target (such as tremor intensity)
    //      This could output a vector for multidimensional targets however other parts of the code will need changing
    //x is your input vector (what your controller produces, and where in the input domain you are sampling)
    //time_iteration is ignored here but is passed for the sake of aligning time points 
    //      (could be key for synchronsing the controller and data acquisition/real system)
    //generate_noise and noise_sd are for generating a noisy update. This was important for testing but obviosuly is 
    //      useless in data acqusition. (either remove or set to 0)
    // testvariables are how the test code changes the properties of the simulated system. This is again useless for
    //      real data.

    int noise_option = 1;
    float mean = 0;

    std::normal_distribution<> noise_dist(mean, noise_sd);
    //this is the one chose for testing
    float next_value = testvariables(0)*((testvariables(1)*sin(x(0)+testvariables(2)))+((1-testvariables(1))*sin((2*x(0))+testvariables(2)))) ; //phase
    next_value = next_value * (1+testvariables(3)-(testvariables(4)*(x(1)-testvariables(5))*(x(1)-testvariables(5))));//frequency
    next_value = next_value * (1+testvariables(6)-(testvariables(7)*(x(2)-testvariables(8))*(x(2)-testvariables(8))));//pulse width
    //noise or not?
    if (noise_option==1) {
        next_value = next_value + noise_dist(generate_noise);
    }
    return next_value;
}


///fill the retained test data. Uses equally spaced tests.
test_data initialise_data(int length_dataset, float x1_min, float x1_max, float x2_min, float x2_max, float x3_min, float x3_max, float noise_sd, VectorXd test_variables){

    test_data dataset;
    dataset.time = Eigen::VectorXd::Zero(length_dataset);
    dataset.inputs_samples = Eigen::MatrixXd::Zero(length_dataset, 3);
    dataset.output_samples = Eigen::VectorXd::Zero(length_dataset);
    ///domain
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist1(x1_min, x1_max);
    std::uniform_real_distribution<> dist2(x2_min, x2_max);
    std::uniform_real_distribution<> dist3(x3_min, x3_max);
    for (int i = 0; i < length_dataset; i++){
       dataset.inputs_samples(i,0) = dist1(gen);
       dataset.inputs_samples(i,1) = dist2(gen);
       dataset.inputs_samples(i,2) = dist3(gen);
    }
    ///get_value
    for (int i = 0; i < length_dataset; i++){

        dataset.output_samples(i) = get_result(dataset.inputs_samples(i,all), 0, gen, noise_sd, test_variables);
    }
    return dataset;
}

///updates the time in the test data
test_data uniform_time_increase(test_data dataset, float time_increase){
    for (int i = 0; i < dataset.time.size(); i++){
       dataset.time(i) += time_increase;
    }
    return dataset;
}

gaussian_process build_model(test_data dataset, int model_resolution_1, int model_resolution_2, int model_resolution_3,
                              float x1_min, float x1_max, float x2_min, float x2_max, float x3_min, float x3_max,
                              float model_decay_const){
    // this is the main part of the controller
    // inputs are: 
    // data (exisiting samples)
    // Model domain details (resultion, range)
    // and a decay constant which discounts old samples -> this should be matched to your system


    gaussian_process model;
    //hyper parameters (see kernel details about these. Will have to match to your system)
    float noise = 1;
    float kernel_hyp_l[3] = {0.6,25,75};//{6.3*0.1,104*0.5,250*0.5};
    float kernel_hyp_s[3] = {3,3,3};//{0.2,0.025,0.025};
    ///generate test domain for your model
    int model_pop = model_resolution_1*model_resolution_2*model_resolution_3;
    model.domain = Eigen::MatrixXd::Zero(model_pop,3);
    int counter = 0;
    for (int a = 0; a < model_resolution_1; a++){
    for (int b = 0; b < model_resolution_2; b++){
    for (int c = 0; c < model_resolution_3; c++){
        model.domain(counter,0) = x1_min + (a*(x1_max-x1_min)/(model_resolution_1-1));
        model.domain(counter,1) = x2_min + (b*(x2_max-x2_min)/(model_resolution_2-1));
        model.domain(counter,2) = x3_min + (c*(x3_max-x3_min)/(model_resolution_3-1));
        counter += 1;
    }}}
    model.mew = Eigen::VectorXd::Zero(model_pop);
    model.sd = Eigen::VectorXd::Zero(model_pop);

    ///gp model
    int length_domain = model.domain.rows();
    int length_dataset = dataset.inputs_samples.rows();
    VectorXd kstar(length_dataset);
    float kstarstar;
    VectorXd xstar(3);
    VectorXd x(3);
    MatrixXd bigk(length_dataset,length_dataset);


    ///big K only needs to be calculated once
    for (int i = 0; i < length_dataset; i++){
        for (int j = 0; j < length_dataset; j++){
            if (i == j) { //can probably speed up if remove if and replace with identity matrix but not sure if it would speed things up much
                bigk(i,j) = calc_kernel(dataset.inputs_samples(i,all),dataset.inputs_samples(j,all), kernel_hyp_s, kernel_hyp_l);
                bigk(i,j) = bigk(i,j)+ noise;// * pow(model_decay_const,((dataset.time(i)+(dataset.time(j))/2))) // account for time decay in relevance
            } else {
                bigk(i,j) = calc_kernel(dataset.inputs_samples(i,all),dataset.inputs_samples(j,all), kernel_hyp_s, kernel_hyp_l);
                bigk(i,j) = bigk(i,j); //* pow(model_decay_const,((dataset.time(i)+(dataset.time(j))/2))); // account for time decay in relevance
            }
        }

    }
    //this part can be optimised as inverse is a cost heavy function
    MatrixXd bigk_inv = bigk.inverse();
    ///loop though domain to compare domain entries to data
    for (int i = 0; i < (length_domain); i++){
        xstar = model.domain(i,all);
        for (int j = 0; j < length_dataset; j++){
            x = dataset.inputs_samples(j,all);
            kstar(j) = calc_kernel(xstar, x, kernel_hyp_s, kernel_hyp_l);
            kstar(j) = kstar(j); //* pow(model_decay_const,(dataset.time(j)/2));;// account for time decay in relevance
            }
        kstarstar = calc_kernel(xstar,xstar, kernel_hyp_s, kernel_hyp_l);
        model.mew(i)  = (kstar.transpose() * bigk_inv * dataset.output_samples);
        model.sd(i) = kstarstar - (kstar.transpose()*bigk_inv*kstar);


    }
    return model;
}

///get_next target value
target_stim aquisition_func(gaussian_process model, float aquisition_variable, int model_res_1, int model_res_2, int model_res_3, float random_number){
    // This does the control policy/adaptive optimisation
    // Essentially it minimises a dynamically weighted balance of predicted cost value and uncertainty.
    // this can be re written as a maximisation problem if desired (see min)
    


    ///min -> this could be flipped to a maximisation problem (maximise mew+b*sd or min mew-b*sd)
    VectorXd model_minimum = model.mew - aquisition_variable * model.sd;


    //cout<<"\n"<<model.sd.transpose()<<"\n";
    Eigen::Index   min_index;
    model_minimum.minCoeff(&min_index);
    target_stim next_value; //{model.domain(0,0), model.domain(0,1), model.domain(0,2)};
    //next_value = {model.domain(0,0), model.domain(0,1), model.domain(0,2)}
    VectorXd next_value_before_rand = model.domain(min_index,all);

    VectorXd model_spacing(2);
    model_spacing(0) = abs(model.domain(0,0)-model.domain(last,0))/model_res_1;
    model_spacing(1) = abs(model.domain(0,1)-model.domain(last,1))/model_res_2;//next_value.stim_variable = next_value_before_rand+model_spacing*random_number;
    next_value.stim_variable = next_value_before_rand;
    next_value.expected_sig = model.sd[min_index];
    next_value.expected_mew = model.mew[min_index];
    ///!!!!!should include safety here to clip within designated regions if implemented for anything with saftey concerns!!!!!!!
    return next_value;
}


float evaluate_information(float expected_sig, float expected_mew, float true_change_trem_intensity, float decay_k, float K){
    //This function is used to dynamically change the Beta parameter


    extern float diff_exp;///////////////

    float new_aq_var;

    float diff = true_change_trem_intensity-expected_mew;


    float limited_sig; //places limits on the model uncertainy
    if (expected_sig<0.01){
        limited_sig = 0.01;
    }else if (expected_sig>10){
        limited_sig = 10;
    }else{
        limited_sig = expected_sig;
    }

    //places limits on the measured uncertainty
    if (diff<0.0001){
        diff = 0.0001;
        //cout<<"_trimmeddifflow_";
    }else if (diff>100){
        //cout<<"_trimmeddiffhigh_";
        diff = 100;
    }
    //played alot with this function: (you can pick anything from the pdf and implement it.)
    diff_exp = K*pow((abs(diff)/(limited_sig)),4) + (decay_k*diff_exp);

    //diff_exp =  k*abs(diff); //y-mu

    new_aq_var = diff_exp;

    //cout<<"next beta"<<new_aq_var;
    return new_aq_var;
}



///select which value to replace
// Key part of TV-GP-UCB (remove samples as they are old)
// Here we want to replace old samples (as they are less likely to be relevant)
// But also want to replace values which are similar to the current sample
// (i.e we dont want to accumulate samples in one area of the input domain)
// atm this is weighted 1:1 (when selection is calculated)
// old samples are already discounted in selecting the next sample point 
int replace_next_value(test_data dataset, gaussian_process model, target_stim aquisition_target, float change_in_tremor_intensity_){
    int n=0;
    //float current_max_time = dataset.time(0);
    //similarity with new value = smaller so inverse
    float dot = abs(dataset.inputs_samples(0,0)-aquisition_target.stim_variable(0));
    dot += abs(dataset.inputs_samples(0,1)-aquisition_target.stim_variable(1));
    dot += abs(dataset.inputs_samples(0,2)-aquisition_target.stim_variable(2));
    float selection = -1*dot - 1*dataset.time(0);
    float max_selection = selection;
    for (int i = 1; i < dataset.time.size(); i++){
        dot = abs(dataset.inputs_samples(i,0)-aquisition_target.stim_variable(0));
        dot += abs(dataset.inputs_samples(i,1)-aquisition_target.stim_variable(1));
        dot += abs(dataset.inputs_samples(i,2)-aquisition_target.stim_variable(2));
        selection = -1*dot - 1*dataset.time(i);
        if (selection > max_selection){
            n = i;
            max_selection = selection;
        }
    }


    //pick replacement

    return n;
}





