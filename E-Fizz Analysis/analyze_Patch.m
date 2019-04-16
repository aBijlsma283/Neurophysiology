function analyze_Patch ()
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_Patch is a mother function that is used to analyze VC and CC
% patch clamp data. The function works for HEKA patch clamp setups EPC 9
% and EPC 10.

% analyze_patch calls other analyze_* functions depending on input. This
% version has the possibility to analyze the ST100 and VCstep protocols.
% Next version will also support the ST10, ST50, ST200 and CC_C functions.

% Protocol command entry:
% VC-Sawtooth (200ms, 100ms, 50ms, 10ms) = 'VCsaw'
% VC-step (all step protocols) = 'VCstep'
% VC-Activate/Inactivate (single analysis only) = 'VCai'

tic
warning off

close all;

addpath(genpath('.'))

%% Protocol extraction

answer = input('Which protocol is used? (VCsaw/VCstep/VCai): ','s');

if strcmp(answer, 'VCsaw')
    protocol = 'VCsaw';
elseif strcmp(answer, 'VCstep')
    protocol = 'VCstep';
elseif strcmp(answer, 'VCai')
    protocol = 'VCai';
else
    disp('Unknown protocol... Analyze_Patch only support VCsaw/VCstep/VCai')
    return
end

clearvars answer
    
%% Visual analysis
% plots every experiments so 'bad' experiments can be filtered out first.

answer = input('Do you want to visual analyse the experiments first? (Yes/No): ','s');

if strcmp(answer, 'Yes') || strcmp(answer, 'yes')
    visualfiles = dir ('*.mat'); % Retrieves all the .mat files data in the current directory.
    for ctr = 1:size(visualfiles,1)
        analyze_Visual(visualfiles(ctr).name(1:end-4))
    end
end

clearvars visualfiles answer

%% Single data analysis
% This part will loop over all files and performs the analysis by the protocol specific
% function (analyze_VCsaw or analyze_VCstep). More information of these
% function can be found by putting 'help analyze_VCsaw' in your command
% window. Each loop will result in an analyzed .mat file, .fig file and the
% original data

answer = input('Do you want to perform single experiment analysis? (Yes/No): ','s');

if strcmp(answer, 'Yes') || strcmp(answer, 'yes')
    foldername = input('Provide the foldername in which the analyzed files will be placed: ','s');
    files = dir('*.mat');
    for ctr = 1:size(files,1)
        name = strsplit(files(ctr).name,'_');
        if strcmp(protocol,'VCsaw') == 1
            try
                analyze_VCsaw(files(ctr).name(1:end-4), foldername)
            catch
                disp(['ERROR: ' files(ctr).name(1:end-4)])
                if exist('Errors','dir') == 7
                    movefile([files(ctr).name(1:end-4) '.mat'] ,'Errors')
                else
                    mkdir('Errors');
                    movefile([files(ctr).name(1:end-4) '.mat'] ,'Errors')
                end
            end
        elseif strcmp(protocol,'VCstep') == 1
            try
                analyze_VCstep(files(ctr).name(1:end-4), foldername)
            catch
                disp(['ERROR: ' files(ctr).name(1:end-4)])
                if exist('Errors','dir') == 7
                    movefile([files(ctr).name(1:end-4) '.mat'] ,'Errors')
                else
                    mkdir('Errors');
                    movefile([files(ctr).name(1:end-4) '.mat'] ,'Errors')
                end
            end
        elseif strcmp(protocol,'VCai') == 1
            try
                analyze_VCai(files(ctr).name(1:end-4), foldername)
            catch
                disp(['ERROR: ' files(ctr).name(1:end-4)])
                if exist('AnalysisErrors','dir') == 7
                    movefile([files(ctr).name(1:end-4) '.mat'] ,'Errors')
                else
                    mkdir('AnalysisErrors');
                    movefile([files(ctr).name(1:end-4) '.mat'] ,'Errors')
                end
            end
        end
    end
    
    cd(foldername)
    analyzedfiles = dir('*.mat');
    disp('Plotting single figures.....')
    
    for ctr = 1:size(analyzedfiles,1)
        try
            analyze_Plotter(analyzedfiles(ctr).name(1:end-4), protocol, foldername)
        catch
            disp(['ERROR: ' analyzedfiles(ctr).name(1:end-4)])
            if exist('PlotErrors','dir') == 7
                movefile([analyzedfiles(ctr).name(1:end-4) '.mat'] ,'PlotErrors')
            else
                mkdir('PlotErrors');
                movefile([analyzedfiles(ctr).name(1:end-4) '.mat'] ,'PlotErrors')
            end
        end
    end
end

clearvars answer

%% Group analysis 
% Retrieves all the analyzed .mat files from the savefolder so it can do a
% groupanalysis. 

answer1 = input('Do you want to perform group analysis? (Yes/No): ','s');
if strcmp(answer1, 'Yes') || strcmp(answer1, 'yes')
    answer2 = input('Do you want to compare two groups with each other (Yes/No): ','s');
    if strcmp(answer2, 'No') || strcmp(answer2, 'no')
        groupname = input('Provide a name for this analysis: ','s');
        if strcmp(protocol,'VCsaw') == 1
            try
                disp('VCsaw group analysis initiated')
                analyze_VCsaw_group(groupname)
            catch
                disp('ERROR: VCsaw group analysis failed...')
            end
        elseif strcmp(protocol,'VCstep') == 1
            try
                disp('VCstep group analysis initiated')
                analyze_VCstep_group(groupname)
            catch
                disp('ERROR: VCstep group analysis failed...')
            end
        elseif strcmp(protocol,'VCai') == 1
            try
                disp('VCai group analysis initiated')
                analyze_VCai_group(groupname)
            catch
                disp('ERROR: VCai group analysis failed...')
            end
        end
    elseif strcmp(answer2, 'Yes') || strcmp(answer2, 'yes')
        groupname = input('Provide an experimentname: ','s');
        group1 = input('Group 1: ','s');
        group2 = input('Group 2: ','s');
        mode = input('Are the two groups paired? (Yes/No): ','s');
        if strcmp(protocol,'VCsaw') == 1
            try
                disp('VCsaw group analysis initiated')
                analyze_VCsaw_group(groupname,group1,group2,mode)
            catch
                disp('ERROR: VCsaw group analysis failed...')
            end
        elseif strcmp(protocol,'VCstep') == 1
            try
                disp('VCstep group analysis initiated')
                analyze_VCstep_group(groupname,group1,group2,mode)
            catch
                disp('ERROR: VCstep group analysis failed...')
            end
        elseif strcmp(protocol,'VCai') == 1
            try
                disp('VCstep group analysis initiated')
                analyze_VCai_group(groupname,group1,group2,mode)
            catch
                disp('ERROR: VCai group analysis failed...')
            end
        end
    end
end
toc