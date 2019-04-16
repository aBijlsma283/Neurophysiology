function analyze_VCai_group(filename,group1,group2)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_VCstep_group is a function used for a voltage clamp sawtooth protocol
% analysis; (magneto version)

% This function retrives all the files pretagged with "analyzed_" and
% performs group analysis on the variables calculated by the analyze_ST functions. 

close all;

%% import the data from the current folder 
if nargin == 1
    filesgroup1 = dir ('analyzed_VCai*.mat');
elseif nargin == 4
    filesgroup1 = dir (['analyzed_VCai*' '_' (group1) '*.mat']);
    filesgroup2 = dir (['analyzed_VCai*' '_' (group2) '*.mat']);
end

%% Gathering of data

if nargin == 1
    for ctr = 1:size(filesgroup1,1)
        load (filesgroup1(ctr).name, 'output')
        gdata.step.preprom(:,1) = output.data.plots.step.pks(:,1);
        gdata.step.preprom(gdata.step.preprom(:,1) == 0) = NaN;
        gdata.step.prom(ctr,1) = nanmean(gdata.step.preprom(:,1));
        gdata.constant.preprom(:,1) = output.data.plots.constant.pks(:,1);
        gdata.constant.preprom(gdata.constant.preprom(:,1) == 0) = NaN;
        gdata.constant.prom(ctr,1) = nanmean(gdata.constant.preprom(:,1));
       
        gdata.step.IVgroup1(ctr,:) = gdata.step.preprom(:,1);
        gdata.step.stimgroup1(ctr,:) = round(output.Vmon.raw(:,1025));
        gdata.constant.IVgroup1(ctr,:) = gdata.constant.preprom(:,1);
        gdata.constant.stimgroup1(ctr,:) = round(output.Vmon.raw(:,1025));
                
        output.data.plots.step.pksindex(isnan(output.data.plots.step.pksindex(:,1)),1) = 1001;
        gdata.step.pretiming(:,1) = output.data.time(output.data.plots.step.pksindex(:,1)-1000);
        gdata.step.pretiming(gdata.step.pretiming(:,1) == 0) = NaN;
        gdata.step.timing(ctr,1) = nanmean(gdata.step.pretiming(:,1));
        output.data.plots.constant.pksindex(isnan(output.data.plots.constant.pksindex(:,1)),1) = 2001;
        output.data.plots.constant.pksindex(output.data.plots.constant.pksindex(:,1) < 2000) = 2001;
        gdata.constant.pretiming(:,1) = output.data.time(output.data.plots.constant.pksindex(:,1)-2000);
        gdata.constant.pretiming(gdata.constant.pretiming(:,1) == 0) = NaN;
        gdata.constant.timing(ctr,1) = nanmean(gdata.constant.pretiming(:,1));
        
        gdata.step.events(ctr,1) = mean(output.data.plots.step.events(:,1));
        gdata.constant.events(ctr,1) = mean(output.data.plots.constant.events(:,1));
        
        gdata.step.preperdiffgroup1(:,1) = output.data.plots.step.pksindex(:,2)-output.data.plots.step.pksindex(:,1);
        gdata.step.preperdiffgroup1(gdata.step.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.step.preperdiffgroup1(:,1)))) == 0
            gdata.step.perdiffgroup1(ctr,1) = output.data.time(round(nanmean(gdata.step.preperdiffgroup1(:,1))));
        else
            gdata.step.perdiffgroup1(ctr,1) = NaN;
        end
        gdata.constant.preperdiffgroup1(:,1) = output.data.plots.constant.pksindex(:,2)-output.data.plots.constant.pksindex(:,1);
        gdata.constant.preperdiffgroup1(gdata.constant.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.constant.preperdiffgroup1(:,1)))) == 0
            gdata.constant.perdiffgroup1(ctr,1) = output.data.time(round(nanmean(gdata.constant.preperdiffgroup1(:,1))));
        else
            gdata.constant.perdiffgroup1(ctr,1) = NaN;
        end
    end
    
elseif nargin == 3
    
    for ctr = 1:size(filesgroup1,1)
        load (filesgroup1(ctr).name, 'output')
        gdata.step.preprom(:,1) = output.data.plots.step.pks(:,1);
        gdata.step.preprom(gdata.step.preprom(:,1) == 0) = NaN;
        gdata.step.prom(ctr,1) = nanmean(gdata.step.preprom(:,1));
        gdata.constant.preprom(:,1) = output.data.plots.constant.pks(:,1);
        gdata.constant.preprom(gdata.constant.preprom(:,1) == 0) = NaN;
        gdata.constant.prom(ctr,1) = nanmean(gdata.constant.preprom(:,1));
       
        gdata.step.IVgroup1(ctr,:) = gdata.step.preprom(:,1);
        gdata.step.stimgroup1(ctr,:) = round(output.Vmon.raw(:,1025));
        gdata.constant.IVgroup1(ctr,:) = gdata.constant.preprom(:,1);
        gdata.constant.stimgroup1(ctr,:) = round(output.Vmon.raw(:,1025));
                
        output.data.plots.step.pksindex(isnan(output.data.plots.step.pksindex(:,1)),1) = 1001;
        gdata.step.pretiming(:,1) = output.data.time(output.data.plots.step.pksindex(:,1)-1000);
        gdata.step.pretiming(gdata.step.pretiming(:,1) == 0) = NaN;
        gdata.step.timing(ctr,1) = nanmean(gdata.step.pretiming(:,1));
        output.data.plots.constant.pksindex(isnan(output.data.plots.constant.pksindex(:,1)),1) = 2001;
        output.data.plots.constant.pksindex(output.data.plots.constant.pksindex(:,1) < 2000) = 2001;
        gdata.constant.pretiming(:,1) = output.data.time(output.data.plots.constant.pksindex(:,1)-2000);
        gdata.constant.pretiming(gdata.constant.pretiming(:,1) == 0) = NaN;
        gdata.constant.timing(ctr,1) = nanmean(gdata.constant.pretiming(:,1));
        
        gdata.step.events(ctr,1) = mean(output.data.plots.step.events(:,1));
        gdata.constant.events(ctr,1) = mean(output.data.plots.constant.events(:,1));
        
        gdata.step.preperdiffgroup1(:,1) = output.data.plots.step.pksindex(:,2)-output.data.plots.step.pksindex(:,1);
        gdata.step.preperdiffgroup1(gdata.step.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.step.preperdiffgroup1(:,1)))) == 0
            gdata.step.perdiffgroup1(ctr,1) = output.data.time(round(nanmean(gdata.step.preperdiffgroup1(:,1))));
        else
            gdata.step.perdiffgroup1(ctr,1) = NaN;
        end
        gdata.constant.preperdiffgroup1(:,1) = output.data.plots.constant.pksindex(:,2)-output.data.plots.constant.pksindex(:,1);
        gdata.constant.preperdiffgroup1(gdata.constant.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.constant.preperdiffgroup1(:,1)))) == 0
            gdata.constant.perdiffgroup1(ctr,1) = output.data.time(round(nanmean(gdata.constant.preperdiffgroup1(:,1))));
        else
            gdata.constant.perdiffgroup1(ctr,1) = NaN;
        end
    end
    
    for ctr = 1:size(filesgroup2,1)
        load (filesgroup2(ctr).name, 'output')
        gdata.step.preprom(:,1) = output.data.plots.step.pks(:,1);
        gdata.step.preprom(gdata.step.preprom(:,1) == 0) = NaN;
        gdata.step.prom(ctr,2) = nanmean(gdata.step.preprom(:,1));
        gdata.constant.preprom(:,1) = output.data.plots.constant.pks(:,1);
        gdata.constant.preprom(gdata.constant.preprom(:,1) == 0) = NaN;
        gdata.constant.prom(ctr,2) = nanmean(gdata.constant.preprom(:,1));
       
        gdata.step.IVgroup2(ctr,:) = gdata.step.preprom(:,1);
        gdata.step.stimgroup2(ctr,:) = round(output.Vmon.raw(:,1025));
        gdata.constant.IVgroup2(ctr,:) = gdata.constant.preprom(:,1);
        gdata.constant.stimgroup2(ctr,:) = round(output.Vmon.raw(:,1025));
                
        output.data.plots.step.pksindex(isnan(output.data.plots.step.pksindex(:,1)),1) = 1001;
        gdata.step.pretiming(:,1) = output.data.time(output.data.plots.step.pksindex(:,1)-1000);
        gdata.step.pretiming(gdata.step.pretiming(:,1) == 0) = NaN;
        gdata.step.timing(ctr,2) = nanmean(gdata.step.pretiming(:,1));
        output.data.plots.constant.pksindex(isnan(output.data.plots.constant.pksindex(:,1)),1) = 2001;
        output.data.plots.constant.pksindex(output.data.plots.constant.pksindex(:,1) < 2000) = 2001;
        gdata.constant.pretiming(:,1) = output.data.time(output.data.plots.constant.pksindex(:,1)-2000);
        gdata.constant.pretiming(gdata.constant.pretiming(:,1) == 0) = NaN;
        gdata.constant.timing(ctr,2) = nanmean(gdata.constant.pretiming(:,1));
        
        gdata.step.events(ctr,2) = mean(output.data.plots.step.events(:,1));
        gdata.constant.events(ctr,2) = mean(output.data.plots.constant.events(:,1));
        
        gdata.step.preperdiffgroup1(:,1) = output.data.plots.step.pksindex(:,2)-output.data.plots.step.pksindex(:,1);
        gdata.step.preperdiffgroup1(gdata.step.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.step.preperdiffgroup1(:,1)))) == 0
            gdata.step.perdiffgroup1(ctr,2) = output.data.time(round(nanmean(gdata.step.preperdiffgroup1(:,1))));
        else
            gdata.step.perdiffgroup1(ctr,2) = NaN;
        end
        gdata.constant.preperdiffgroup1(:,1) = output.data.plots.constant.pksindex(:,2)-output.data.plots.constant.pksindex(:,1);
        gdata.constant.preperdiffgroup1(gdata.constant.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.constant.preperdiffgroup1(:,1)))) == 0
            gdata.constant.perdiffgroup1(ctr,2) = output.data.time(round(nanmean(gdata.constant.preperdiffgroup1(:,1))));
        else
            gdata.constant.perdiffgroup1(ctr,2) = NaN;
        end
    end
end
keyboard
%% boxplots (figure 1)

save(['completegroup_VCai_' filename], 'gdata', 'output')