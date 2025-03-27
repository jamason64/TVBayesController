#include <iostream>
#include "gplibrary.hpp" // contains the Gaussian process functions
#include "test_functions.hpp" // contains the code for generating test functions


using namespace std;
using Eigen::MatrixXd; // set up for Eigen (matrix&vectors)
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::last;

///for csv > records tests and is read by jupyter notebook
#include <fstream>
#include <string>
//retained data going to csv
string datatime; //age of the data sample in cycle number
string dataphase; //sample data dimension 1
string datavoltage; // sample data dimension 2 (became frequency
string datapulseduration;//sample data dimension 3
string dataoutput;  //output value of sample data
//model variable going to csv
string modelphase; //dimension 1
string modelvoltage; //dimension 2
string modelpulseduration; //dimension 3
string modelmew; //model's mean output
string modelsd; //model's standard deviation output
///






const int remember_n_stims =7; //The number of previous samples that feed the GP model

float number_of_optimisation_iterations = 50; //for testing, number of sample,model cycles

const int model_resolution_phase =  6; // dimensions for the model's parameter space
const int model_resolution_freq = 3;
const int model_resolution_pulse = 6;


const float phase_rad_min = 0; //max,min for parameter space dimensions
const float phase_rad_max = 6.28;
const float freq_hz_min = 52;
const float freq_hz_max = 156;
const float pulse_duration_min_us = 150;
const float pulse_duration_max_us = 400;





//GP model/controller constants
float aquisition_variable = 0.009; //
float diff_exp = 0; //
float aquisition_proportional_const = 1; //
float decay_constant = 0.45; //

const float standard_deviation_noise = 0.1;
const float model_time_decay = 0.1;

int replacement_index;



//these are for controlling the testing
int test_choice = 1; // 1 = step, 3 = random drift
float test_step_size = 0.000001; //0 to 1 division of given ranges fir each step (gives rate of drift)
VectorXd test_var;//before step []
ranges var_range;




///for lib uses///
test_data retained_dataset;
gaussian_process model;
target_stim next_stimulation;
float change_in_tremor_intensity;


int main() {
    //set up repeatable random generation
    std::random_device rd;
    srand((15));//***** This number is what is changed to get different test trials. (I kept development and test seeds seperate)
    int random = rand(); // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(random);
    std::uniform_real_distribution<> distnv(-0.5, 0.5);




    ///for the random generation of the objective function, these functions&variables can be found in the write=-up
    var_range.phase_m[0] = 0;
    var_range.phase_m[1] = 1;
    var_range.phase_a[0] = 0;
    var_range.phase_a[1] = 1;
    var_range.phase_b[0] = 0;
    var_range.phase_b[1] = 6.28;

    var_range.freq_a[0] = 0;
    var_range.freq_a[1] = 0.2;
    var_range.freq_b[0] = 0.000005;
    var_range.freq_b[1] = 0.00005;
    var_range.freq_s[0] = 50;
    var_range.freq_s[1] = 200;

    var_range.pw_a[0] = 0;
    var_range.pw_a[1] = 0.2;
    var_range.pw_b[0] = 0.000001;
    var_range.pw_b[1] = 0.00001;
    var_range.pw_s[0] = 150;
    var_range.pw_s[1] = 400;

    test_var = Eigen::VectorXd::Zero(9);
    test_var = get_test_var(test_var, var_range, 1, gen, 1); // generate sample 1's objective function
    //cout<<"\n step \n"<<test_var<<" \n";
    ///end of test stuff






    std::ofstream myfile;//csv line


    ///initialisation of remembered samples feeding model (equivalent of in-clinic set-up)
    retained_dataset = initialise_data(remember_n_stims,
                                        phase_rad_min, phase_rad_max,
                                        freq_hz_min, freq_hz_max,
                                         pulse_duration_min_us, pulse_duration_max_us,
                                         standard_deviation_noise, test_var);



    //The csv lines feed the jupyter notebook which allows for interpretation of tests.
    ///csv
    myfile.open ("data.csv");
    myfile<< "datatime,dataphase,datafrequency,datapulseduration,dataoutput,";
    myfile<<  "modelphase,modelfrequency,modelpulseduration,modelmew,modelsd,beta,";
    myfile<<  "phase_m,phase_a,phase_b,freq_a,freq_b,freq_s,pw_a,pw_b,pw_s,\n";
    ///csv














    ///iterations
    //each cycle is a sample
    for (int run_n = 0; run_n < number_of_optimisation_iterations; run_n++){

        switch (test_choice){
            case(1):{
                if (run_n == 25){test_var = get_test_var(test_var, var_range, 1, gen, 1); cout<<"\n step \n"<<test_var<<" \n";}
            break;}
            case(2):{
                test_var = get_test_var(test_var, var_range, 2, gen, test_step_size);
            break;}
        }


        model = build_model(retained_dataset, model_resolution_phase, model_resolution_freq, model_resolution_pulse,
                              phase_rad_min, phase_rad_max,
                               freq_hz_min, freq_hz_max,
                                pulse_duration_min_us, pulse_duration_max_us,
                                model_time_decay);
        //cout<<"\n"<<model.mew;


        next_stimulation = aquisition_func(model, aquisition_variable, model_resolution_phase, model_resolution_freq, model_resolution_pulse, distnv(gen));
        //cout<<"\n" <<"balance m"<< next_stimulation.expected_mew <<"balance s "<< next_stimulation.expected_sig<<"balance beta "<<aquisition_variable;

        ///just to get the stationary results:
        //next_stimulation.stim_variable(0) = 5.25;
        //next_stimulation.stim_variable(1) = 52;
        //next_stimulation.stim_variable(2) = 250;

        change_in_tremor_intensity = get_result(next_stimulation.stim_variable, run_n, gen, standard_deviation_noise, test_var);
        //cout<<"x"<<change_in_tremor_intensity;


        aquisition_variable = evaluate_information(next_stimulation.expected_sig, next_stimulation.expected_mew, change_in_tremor_intensity, decay_constant, aquisition_proportional_const);
        //0.01; //

        //cout<<"diff_exp"<<diff_exp<<"\n";
        //cout << "beta " <<aquisition_variable<< "\n";
        //cout << "s " <<next_stimulation.expected_sig<< "\n";
        //cout << "diff " <<abs(next_stimulation.expected_mew-change_in_tremor_intensity)<< "\n";
        //cout<<"iti \n";
        replacement_index = replace_next_value(retained_dataset, model,
                                                next_stimulation, change_in_tremor_intensity); ///to do ***** use this///
        retained_dataset = uniform_time_increase(retained_dataset, 1);
        //cout<<"ri \n";
        retained_dataset.inputs_samples(replacement_index,0) = next_stimulation.stim_variable[0];
        retained_dataset.inputs_samples(replacement_index,1) = next_stimulation.stim_variable[1];
        retained_dataset.inputs_samples(replacement_index,2) = next_stimulation.stim_variable[2];
        retained_dataset.output_samples(replacement_index) = change_in_tremor_intensity;
        retained_dataset.time(replacement_index) = 0;

        //cout << "it " << run_n <<"\n";


        ///csv

        for (int csvlinen=0; csvlinen<remember_n_stims; csvlinen++){
            myfile << to_string(retained_dataset.time(csvlinen))+",";
            myfile << to_string(retained_dataset.inputs_samples(csvlinen, 0))+",";
            myfile << to_string(retained_dataset.inputs_samples(csvlinen, 1))+",";
            myfile << to_string(retained_dataset.inputs_samples(csvlinen, 2))+",";
            myfile << to_string(retained_dataset.output_samples(csvlinen))+",";
            myfile << to_string(model.domain(csvlinen, 0))+",";
            myfile << to_string(model.domain(csvlinen, 1))+",";
            myfile << to_string(model.domain(csvlinen, 2))+",";
            myfile << to_string(model.mew(csvlinen))+",";
            myfile << to_string(model.sd(csvlinen))+",";
            myfile << to_string(aquisition_variable)+",";
            myfile << to_string(test_var(0))+",";
            myfile << to_string(test_var(1))+",";
            myfile << to_string(test_var(2))+",";
            myfile << to_string(test_var(3))+",";
            myfile << to_string(test_var(4))+",";
            myfile << to_string(test_var(5))+",";
            myfile << to_string(test_var(6))+",";
            myfile << to_string(test_var(7))+",";
            myfile << to_string(test_var(8))+",\n";
        }
        for (int csvlinen=remember_n_stims; csvlinen<(model_resolution_phase*model_resolution_pulse*model_resolution_freq); csvlinen++){
            myfile << "NaN,";
            myfile << "NaN,";
            myfile << "NaN,";
            myfile << "NaN,";
            myfile << "NaN,";
            myfile << to_string(model.domain(csvlinen, 0))+",";
            myfile << to_string(model.domain(csvlinen, 1))+",";
            myfile << to_string(model.domain(csvlinen, 2))+",";
            myfile << to_string(model.mew(csvlinen))+",";
            myfile << to_string(model.sd(csvlinen))+",";
            myfile << to_string(aquisition_variable)+",";
            myfile << to_string(test_var(0))+",";
            myfile << to_string(test_var(1))+",";
            myfile << to_string(test_var(2))+",";
            myfile << to_string(test_var(3))+",";
            myfile << to_string(test_var(4))+",";
            myfile << to_string(test_var(5))+",";
            myfile << to_string(test_var(6))+",";
            myfile << to_string(test_var(7))+",";
            myfile << to_string(test_var(8))+",\n";
        }

        ///csv

    }
    ///csv
    myfile.close();
    ///csv
}
