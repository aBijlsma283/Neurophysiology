function analyze_Visual(filename)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_Visual is part of the mother function 'analyze_Patch' and is a function 
% used for a VC/CC visual analysis.

% This function retrieves the data exported from HEKA patch clamp setups
% EPC 9 and EPC 10 double. It puts the data into matlab and performs data
% analysis. The function plots every experiments so good/bad experiments
% can be seperated.

close all; 
warning off
load ([filename '.mat']);

files = whos ('Trace*'); % lists the current variables, in long format form
output.file.name = filename;
disp (['analyze: ' output.file.name]) % displays what you are doing, the VC_Step_Analysis

%%

for Tc = 1:size(files,1) % loops until every trace is used, amount of traces is equal to the size(files,1). As it lists all the traces.
    
    Ts = strsplit(files(Tc).name,'_'); % Splits the traces at the '_'. For our example it will split in 'Trace', '4', '7', '1', '1'.
    
    % Organizing and Assigning the split traces. str2double is used to convert
    % the elements in the cell array of character to vectors C to double
    
    output.exp.experiment (Tc) = str2double(Ts{2});
    output.exp.dataset (Tc) = str2double(Ts{3});
    output.exp.sweep (Tc) = str2double(Ts{4});
    output.exp.channel (Tc) = str2double(Ts{5});
    
    output.data.tracename {Tc} = files(Tc).name;
    output.data.values {Tc} = eval(files(Tc).name);
    output.data.sweepduration (Tc) =  output.data.values {Tc}(end,1)-output.data.values{Tc}(1,1);             % in seconds
    output.data.samplingrate (Tc) =  round(size(output.data.values {Tc},1)/output.data.sweepduration (Tc));   % in Hz
end

output.data.time{1} = output.data.values{1}(:,1);
%%

if round(max(output.data.values{2}(:,2)),2) == 0.0700 || round(max(output.data.values{2}(:,2)),2) == 0.0800 || round(max(output.data.values{2}(:,2)),2) == 0.1000 || max(output.exp.channel) == 2 || round(min(output.data.values{2}(:,2)),2) == -0.0800 || round(min(output.data.values{2}(:,2)),2) == -0.0700 || round(min(output.data.values{2}(:,2)),2) == -0.0900
    Ichannel = 1;
    Vchannel = 2;
else
    Ichannel = 3;
    Vchannel = 4;
end

%%
Vctr = 1;
Ictr = 1;

for Tc = 1:size(files,1)
    
    if output.exp.channel(Tc) == Vchannel % only picks from column 2, as we assigned 2 = Vchannel
        output.Vmon.raw(Vctr,:) = 1000.*output.data.values{Tc}(:,2); % Data is V, we want it in mV so we multiply with 1000.
        
        Vctr = Vctr + 1;
        
    elseif output.exp.channel(Tc) == Ichannel % only picks from column 1, as we assigned 1 = Vchannel
        output.Imon.raw(Ictr,:) = 10e8.*output.data.values{Tc}(:,2); % Data is A, we want it in nA so we multiply with 10e8.
        
        Ictr = Ictr + 1;
    end
end

plot(output.data.time{1},output.Imon.raw)
title(output.file.name,'FontSize',12)
xlabel('time [ms]','FontSize',12)
ylabel('Current [nA]','FontSize',12)
set(gcf,'Position', get(0, 'Screensize'));
  
grade = input('1 = Include, 2 = Exclude, 3 = O_o questionable...... Grade = ');

if grade == 2
    mkdir('excluded') ; movefile([filename '.mat'],'excluded')
elseif grade == 3
    mkdir('questionable') ; movefile([filename '.mat'],'questionable')
end


