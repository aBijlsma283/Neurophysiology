close all; 
warning off

filesmat = dir('*.mat');

for ctr = 1:size(filesmat,1)
    name = strsplit(filesmat(ctr).name,'_');
    Trace_changer_help([filesmat(ctr).name(1:end-4) '.mat']);
    
    movefile('temp.mat',[filesmat(ctr).name(1:end-4) '.mat'])
end