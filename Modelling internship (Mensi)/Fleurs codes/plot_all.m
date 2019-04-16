all_files = dir('*.mat');

gdata.C_Exc.gamma = nan(5,15);
gdata.C_Exc.MI_i = nan(5,15);
gdata.C_Exc.MSE_i = nan(5,15);
gdata.C_Exc.MI = nan(5,15);
gdata.C_Exc.qon = nan(5,15); 
gdata.C_Exc.qoff = nan(5,15); 
gdata.C_Exc.MSE = nan(5,15);
gdata.C_Exc.FI = nan(5,15);
gdata.D1_Exc.gamma = nan(5,15);
gdata.D1_Exc.MI_i = nan(5,15);
gdata.D1_Exc.MSE_i = nan(5,15);
gdata.D1_Exc.MI = nan(5,15);
gdata.D1_Exc.qon = nan(5,15); 
gdata.D1_Exc.qoff = nan(5,15); 
gdata.D1_Exc.MSE = nan(5,15);
gdata.D1_Exc.FI = nan(5,15);
gdata.C_Inh.gamma = nan(5,15);
gdata.C_Inh.MI_i = nan(5,15);
gdata.C_Inh.MSE_i = nan(5,15);
gdata.C_Inh.MI = nan(5,15);
gdata.C_Inh.qon = nan(5,15); 
gdata.C_Inh.qoff = nan(5,15); 
gdata.C_Inh.MSE = nan(5,15);
gdata.C_Inh.FI = nan(5,15);
gdata.D1_Inh.gamma = nan(5,15);
gdata.D1_Inh.MI_i = nan(5,15);
gdata.D1_Inh.MSE_i = nan(5,15);
gdata.D1_Inh.MI = nan(5,15);
gdata.D1_Inh.qon = nan(5,15); 
gdata.D1_Inh.qoff = nan(5,15); 
gdata.D1_Inh.MSE = nan(5,15);
gdata.D1_Inh.FI = nan(5,15);

N1 = 0;
N2 = 0;
N3 = 0;
N4 = 0;

for ctr = 1:size(all_files,1)
    disp(ctr)
    load(all_files(ctr).name(1:end-4))
    name = all_files(ctr).name(1:end-4);
    namesplits = strsplit(name,'_');
    AltN = str2double(namesplits{4});
    
    if strcmp(Data.settings.condition,'aCSF')
        if Data.settings.tau == 50      % Inh
            if AltN == 1
               N1 = N1 + 1; 
            end
            gdata.C_Inh.gamma(N1,AltN) = Data.settings.gamma;
            gdata.C_Inh.MI_i(N1,AltN) = Data.Analysis{1,1}.MI_i;
            gdata.C_Inh.MSE_i(N1,AltN) = Data.Analysis{1,1}.MSE_i;
            gdata.C_Inh.MI(N1,AltN) = Data.Analysis{1,1}.MI;
            gdata.C_Inh.qon(N1,AltN) = Data.Analysis{1,1}.qon;
            gdata.C_Inh.qoff(N1,AltN) = Data.Analysis{1,1}.qoff;
            gdata.C_Inh.MSE(N1,AltN) = Data.Analysis{1,1}.MSE;
            gdata.C_Inh.FI(N1,AltN) = Data.Analysis{1,1}.FI;
            gdata.C_Inh.xhat_i{N1,AltN} = Data.Analysis{1,1}.xhat_i;
            gdata.C_Inh.xhatspikes{N1,AltN} = Data.Analysis{1,1}.xhatspikes;
        elseif Data.settings.tau == 250 % Exc
            if AltN == 1
               N2 = N2 + 1; 
            end
            gdata.C_Exc.gamma(N2,AltN) = Data.settings.gamma;
            gdata.C_Exc.MI_i(N2,AltN) = Data.Analysis{1,1}.MI_i;
            gdata.C_Exc.MSE_i(N2,AltN) = Data.Analysis{1,1}.MSE_i;
            gdata.C_Exc.MI(N2,AltN) = Data.Analysis{1,1}.MI;
            gdata.C_Exc.qon(N2,AltN) = Data.Analysis{1,1}.qon;
            gdata.C_Exc.qoff(N2,AltN) = Data.Analysis{1,1}.qoff;
            gdata.C_Exc.MSE(N2,AltN) = Data.Analysis{1,1}.MSE;
            gdata.C_Exc.FI(N2,AltN) = Data.Analysis{1,1}.FI;
            gdata.C_Exc.xhat_i{N2,AltN} = Data.Analysis{1,1}.xhat_i;
            gdata.C_Exc.xhatspikes{N2,AltN} = Data.Analysis{1,1}.xhatspikes;
        end
    elseif strcmp(Data.settings.condition,'D1ago')
        if Data.settings.tau == 50      % Inh
            if AltN == 1
               N3 = N3 + 1; 
            end
            gdata.D1_Inh.gamma(N3,AltN) = Data.settings.gamma;
            gdata.D1_Inh.MI_i(N3,AltN) = Data.Analysis{1,1}.MI_i;
            gdata.D1_Inh.MSE_i(N3,AltN) = Data.Analysis{1,1}.MSE_i;
            gdata.D1_Inh.MI(N3,AltN) = Data.Analysis{1,1}.MI;
            gdata.D1_Inh.qon(N3,AltN) = Data.Analysis{1,1}.qon;
            gdata.D1_Inh.qoff(N3,AltN) = Data.Analysis{1,1}.qoff;
            gdata.D1_Inh.MSE(N3,AltN) = Data.Analysis{1,1}.MSE;
            gdata.D1_Inh.FI(N3,AltN) = Data.Analysis{1,1}.FI;
            gdata.D1_Inh.xhat_i{N3,AltN} = Data.Analysis{1,1}.xhat_i;
            gdata.D1_Inh.xhatspikes{N3,AltN} = Data.Analysis{1,1}.xhatspikes;
        elseif Data.settings.tau == 250 % Exc
            if AltN == 1
               N4 = N4 + 1; 
            end
            gdata.D1_Exc.gamma(N4,AltN) = Data.settings.gamma;
            gdata.D1_Exc.MI_i(N4,AltN) = Data.Analysis{1,1}.MI_i;
            gdata.D1_Exc.MSE_i(N4,AltN) = Data.Analysis{1,1}.MSE_i;
            gdata.D1_Exc.MI(N4,AltN) = Data.Analysis{1,1}.MI;
            gdata.D1_Exc.qon(N4,AltN) = Data.Analysis{1,1}.qon;
            gdata.D1_Exc.qoff(N4,AltN) = Data.Analysis{1,1}.qoff;
            gdata.D1_Exc.MSE(N4,AltN) = Data.Analysis{1,1}.MSE;
            gdata.D1_Exc.FI(N4,AltN) = Data.Analysis{1,1}.FI;
            gdata.D1_Exc.xhat_i{N4,AltN} = Data.Analysis{1,1}.xhat_i;
            gdata.D1_Exc.xhatspikes{N4,AltN} = Data.Analysis{1,1}.xhatspikes;
        end
    end
end

keyboard

%% plotting

