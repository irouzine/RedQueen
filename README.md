# RedQueen
MATLAB codes for the paper on vaccination effect on virus evolution (Rouzine 2025):

OneDwaveNew.m       % Monte Carlo simulation of the model (Fig. 2)

all_vaccine_w_Dz_bz.m % semi–analytic calculation of evolution rate and incidence when vaccine effects overlap (Fig. 4A-C)

vaccine_w_Dz_bz.m    % semi–analytic calculation of evolution rate and incidence when only last vaccine has a memory effect (Fig. 4D)

vaccine.m            % frequent vaccination limit with ovelapping effects (Fig 3).  

scaling_genealogy.m  % testing Eq. S24 in S2 Text with Monte Carlo simulation (Figure  S1).  

recomb_train.m       % Called by scaling_genealogy.m. Monte Carlo evolution program to generate binary genome sequences. 

recomb_train2.m       % Called by testS3andS5.m. A version of recomb_train.m.

main_fitting.m       % Fitting b_y and U_b to data on influenza, no vaccine (S2 Text, Figure S3 ). 

call_fitting.m       % Called by main_fitting. Calculates the mean square diff between V and TMRCA and the experimental values. 

testS3andS5.m        % Testing Eq. S26 for the substitution rate in the case of a uniform distribution of s (S2 Text, Figure S2). 
