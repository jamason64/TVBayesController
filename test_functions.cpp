#include "gplibrary.hpp"
#include "test_functions.hpp"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::last;
using Eigen::all;

VectorXd get_test_var(VectorXd previous, ranges dist_range, int testoption, std::mt19937 gen_rd, float stepsize){
    //previous objective variables, valid ranges, testoption 1 = random generation|2 = random step, random seed, stepsize for random step
    VectorXd locations(9);
    switch(testoption) {
        case 1: {///random generation
            std::uniform_real_distribution<> dist_phase_m(dist_range.phase_m[0], dist_range.phase_m[1]);
            std::uniform_real_distribution<> dist_phase_a(dist_range.phase_a[0], dist_range.phase_a[1]);
            std::uniform_real_distribution<> dist_phase_b(dist_range.phase_b[0], dist_range.phase_b[1]);
            std::uniform_real_distribution<> dist_freq_a(dist_range.freq_a[0], dist_range.freq_a[1]);
            std::uniform_real_distribution<> dist_freq_b(dist_range.freq_b[0], dist_range.freq_b[1]);
            std::uniform_real_distribution<> dist_freq_s(dist_range.freq_s[0], dist_range.freq_s[1]);
            std::uniform_real_distribution<> dist_pw_a(dist_range.pw_a[0], dist_range.pw_a[1]);
            std::uniform_real_distribution<> dist_pw_b(dist_range.pw_b[0], dist_range.pw_b[1]);
            std::uniform_real_distribution<> dist_pw_s(dist_range.pw_s[0], dist_range.pw_s[1]);


            locations(0) = dist_phase_m(gen_rd);
            locations(1) = dist_phase_a(gen_rd);
            locations(2) = dist_phase_b(gen_rd);
            locations(3) = dist_freq_a(gen_rd);
            locations(4) = dist_freq_b(gen_rd);
            locations(5) = dist_freq_s(gen_rd);
            locations(6) = dist_pw_a(gen_rd);
            locations(7) = dist_pw_b(gen_rd);
            locations(8) = dist_pw_s(gen_rd);
        break;}

        case 2: {///+-step, legacy not used
            std::normal_distribution<> plusminus(0,stepsize);

            locations(0) = previous(0) + (plusminus(gen_rd) * abs(dist_range.phase_m[0]-dist_range.phase_m[1]));
            locations(1) = previous(1) + (plusminus(gen_rd) * abs(dist_range.phase_a[0]- dist_range.phase_a[1]));
            locations(2) = previous(2) + (plusminus(gen_rd) * abs(dist_range.phase_b[0]- dist_range.phase_b[1]));
            locations(3) = previous(3) + (plusminus(gen_rd) * abs(dist_range.freq_a[0]- dist_range.freq_a[1]));
            locations(4) = previous(4) + (plusminus(gen_rd) * abs(dist_range.freq_b[0]- dist_range.freq_b[1]));
            locations(5) = previous(5) + (plusminus(gen_rd) * abs(dist_range.freq_s[0]- dist_range.freq_s[1]));
            locations(6) = previous(6) + (plusminus(gen_rd) * abs(dist_range.pw_a[0]- dist_range.pw_a[1]));
            locations(7) = previous(7) + (plusminus(gen_rd) * abs(dist_range.pw_b[0]- dist_range.pw_b[1]));
            locations(8) = previous(8) + (plusminus(gen_rd) * abs(dist_range.pw_s[0]- dist_range.pw_s[1]));

            if  (locations(0)<dist_range.phase_m[0]){locations(0) = dist_range.phase_m[0];}
            else if (locations(0)>dist_range.phase_m[1]){locations(0) = dist_range.phase_m[1];}
            if  (locations(1)<dist_range.phase_a[0]){locations(1) = dist_range.phase_a[0];}
            else if (locations(1)>dist_range.phase_a[1]){locations(1) = dist_range.phase_a[1];}
            if  (locations(2)<dist_range.phase_b[0]){locations(2) = dist_range.phase_b[0];}
            else if (locations(2)>dist_range.phase_b[1]){locations(2) = dist_range.phase_b[1];}

            if  (locations(3)<dist_range.freq_a[0]){locations(3) = dist_range.freq_a[0];}
            else if (locations(3)>dist_range.freq_a[1]){locations(3) = dist_range.freq_a[1];}
            if  (locations(4)<dist_range.freq_b[0]){locations(4) = dist_range.freq_b[0];}
            else if (locations(4)>dist_range.freq_b[1]){locations(4) = dist_range.freq_b[1];}
            if  (locations(5)<dist_range.freq_s[0]){locations(5) = dist_range.freq_s[0];}
            else if (locations(5)>dist_range.freq_s[1]){locations(5) = dist_range.freq_s[1];}

            if  (locations(6)<dist_range.pw_a[0]){locations(6) = dist_range.pw_a[0];}
            else if (locations(6)>dist_range.pw_a[1]){locations(6) = dist_range.pw_a[1];}
            if  (locations(7)<dist_range.pw_b[0]){locations(7) = dist_range.pw_b[0];}
            else if (locations(7)>dist_range.pw_b[1]){locations(7) = dist_range.pw_b[1];}
            if  (locations(8)<dist_range.pw_s[0]){locations(8) = dist_range.pw_s[0];}
            else if (locations(8)>dist_range.pw_s[1]){locations(8) = dist_range.pw_s[1];}

        break;}

        case 3: {///+-step (drift test)
            std::uniform_int_distribution<> plusminus(0,1);

            locations(0) = previous(0) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.phase_m[0]-dist_range.phase_m[1]) * stepsize);
            locations(1) = previous(1) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.phase_a[0]- dist_range.phase_a[1]) * stepsize);
            locations(2) = previous(2) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.phase_b[0]- dist_range.phase_b[1]) * stepsize);
            locations(3) = previous(3) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.freq_a[0]- dist_range.freq_a[1]) * stepsize);
            locations(4) = previous(4) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.freq_b[0]- dist_range.freq_b[1]) * stepsize);
            locations(5) = previous(5) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.freq_s[0]- dist_range.freq_s[1]) * stepsize);
            locations(6) = previous(6) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.pw_a[0]- dist_range.pw_a[1]) * stepsize);
            locations(7) = previous(7) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.pw_b[0]- dist_range.pw_b[1]) * stepsize);
            locations(8) = previous(8) + (pow(-1,plusminus(gen_rd)) * abs(dist_range.pw_s[0]- dist_range.pw_s[1]) * stepsize);

            if  (locations(0)<dist_range.phase_m[0]){locations(0) = dist_range.phase_m[0];}
            else if (locations(0)>dist_range.phase_m[1]){locations(0) = dist_range.phase_m[1];}
            if  (locations(1)<dist_range.phase_a[0]){locations(1) = dist_range.phase_a[0];}
            else if (locations(1)>dist_range.phase_a[1]){locations(1) = dist_range.phase_a[1];}
            if  (locations(2)<dist_range.phase_b[0]){locations(2) = dist_range.phase_b[0];}
            else if (locations(2)>dist_range.phase_b[1]){locations(2) = dist_range.phase_b[1];}

            if  (locations(3)<dist_range.freq_a[0]){locations(3) = dist_range.freq_a[0];}
            else if (locations(3)>dist_range.freq_a[1]){locations(3) = dist_range.freq_a[1];}
            if  (locations(4)<dist_range.freq_b[0]){locations(4) = dist_range.freq_b[0];}
            else if (locations(4)>dist_range.freq_b[1]){locations(4) = dist_range.freq_b[1];}
            if  (locations(5)<dist_range.freq_s[0]){locations(5) = dist_range.freq_s[0];}
            else if (locations(5)>dist_range.freq_s[1]){locations(5) = dist_range.freq_s[1];}

            if  (locations(6)<dist_range.pw_a[0]){locations(6) = dist_range.pw_a[0];}
            else if (locations(6)>dist_range.pw_a[1]){locations(6) = dist_range.pw_a[1];}
            if  (locations(7)<dist_range.pw_b[0]){locations(7) = dist_range.pw_b[0];}
            else if (locations(7)>dist_range.pw_b[1]){locations(7) = dist_range.pw_b[1];}
            if  (locations(8)<dist_range.pw_s[0]){locations(8) = dist_range.pw_s[0];}
            else if (locations(8)>dist_range.pw_s[1]){locations(8) = dist_range.pw_s[1];}

        break;}

}







    return locations;
};
