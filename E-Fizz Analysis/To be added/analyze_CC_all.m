function analyze_CC_all
% Analyze_CC_all is a function in the intracellular e-fizz toolbox. 
%
% The function retrieves all the .mat data in the current directory. It 
% assumes that the user exported them from HEKA into MATLAB. The function
% then loops over all the data to perform data analysis using analyze_cc
% function. The assumption here is that the data is coming from current
% clamp experiments. For more details on the analyze_CC please enter in
% your command window:
% help analyze_CC
%
% Sample entry:
% analyze_CC_all

% This code uses a function from signal processing toolbox. 

% Tansu Celikel -- celikel@neurophysiology.nl
% V1. 16.06.2016


allexp = dir (['*mat']);

for lp = 1:size(allexp,1)

    disp(['Working on ' allexp(lp).name ' >>> Please wait! <<<'])
    analyze_CC([allexp(lp).name(1:end-4)])

end
