README.txt
==========
Bayesian Optimisation as a Controller for Peripheral Stimulation
James Mason 2023
University of Oxford


Project Overview
----------------
This repository contains C++ code and accompanying files for a time-varying Bayesian optimisation controller. 
It is originally intended to help select stimulation parameters for minimising tremor using a neuromodulation device. 


This software is provided under the MIT Licence.
In addition, the author kindly requests that any commercial use be disclosed to the author (James Mason).
If you use this code in academic work, please cite the GitHub repository: https://github.com/jamason64/TVBayesController.git

Files and Directories
---------------------
1. gplibrary.hpp / gplibrary.cpp
   - The library/controller

2. toplevel.cpp
   - A usage example of the library
   - Set-up for testing as in write-up

3. test_functions.hpp / test_functions.cpp
   - Provide simulated objective functions for testing

4. interpretcsv.ipynb
   - Jupyter notebook for analysing CSV outputs of toplevel.cpp

5. data.csv
   - output of toplevel.cpp.

6. James_Mason_4YP_Submission_160523.pdf
   - Report explaining theory and testing results.

7. README.txt
   - This file.

Functionality
-------------
- Finds a optimal outcome for a changing black box system by continually refining a model mapping stimulation parameters to outcomes (e.g. tremor reduction).
Key features over standard bayesian optimisation: 
- Adjusts exploration bias with the difference between expected and true information gain/uncertainty
- Uses a "sliding window" or "forgetting" approach so that older data is less dominant.

How to Use
----------
Your best bet is to ask chatgpt for how to get this code running on your system. 
toplevel.cpp should run as is integrated with the gplibrary code and the test_functions. 
See toplevel.cpp and the minimal implementation in gplibrary.cpp for how to implment with a real system.



Further Details
---------------
- See James_Mason_4YP_Submission_160523.pdf for the full background and analysis.
- Feel free to contact me at james.w.mason@gmail.com for any more information.

End of README
