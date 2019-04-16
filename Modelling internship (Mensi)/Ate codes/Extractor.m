%% 

CInhfiles = dir('*C_Inh_normal.mat');
CExcfiles = dir('*C_Exc_normal.mat');
D1Inhfiles = dir('*D1_Inh_normal.mat');
D1Excfiles = dir('*D1_Exc_normal.mat');

for ctrI = 1:size(CInhfiles,1)
    load(CInhfiles(ctrI).name(1:end-4))
    gdata.Cinh.mspikestimes{ctrI} = model_params.m_spiketimes*20;
    gdata.Cinh.rspikestimes{ctrI} = model_params.r_spiketimes*20;
%     for i = 1:size(Data.Analysis,2)
%         gdata.Cinh.Mi(ctrI,i) = Data.Analysis{1,i}.MI_i;
%         gdata.Cinh.Fr(ctrI,i) = Data.firing_rate(1,i);
%     end
%     gdata.Cinh.param(ctrI,1:7) = model_params.param;
%     gdata.Cinh.events(ctrI,1:4) = model_params.r_spikevents(2,1:4);
%     gdata.Cinh.spikecount(ctrI) =  model_params.r_spikecount;
%     gdata.Cinh.eta(ctrI,:) =  model_params.eta;
%     gdata.Cinh.gamma(ctrI,:) =  model_params.gamma;
%     gdata.Cinh.spikeshape(ctrI,:) =  model_params.spikeshape;
%     gdata.Cinh.G2(ctrI) = model_params.G2;
%     gdata.Cinh.G5(ctrI) = model_params.G5;
    gdata.Cinh.G10(ctrI) = model_params.G10;
%     gdata.Cinh.Mc(ctrI) = model_params.Mc;
    clearvars -except ctrI CInhfiles CExcfiles gdata D1Inhfiles D1Excfiles
end
% gdata.Cinh.param(:,6) = gdata.Cinh.param(:,6) + gdata.Cinh.Mc';

disp('CInh done')
for ctrE = 1:size(CExcfiles,1)
    load(CExcfiles(ctrE).name(1:end-4))
    gdata.Cexc.mspikestimes{ctrE} = model_params.m_spiketimes*20;
    gdata.Cexc.rspikestimes{ctrE} = model_params.r_spiketimes*20;
%     for i = 1:size(Data.Analysis,2)
%         gdata.CExc.Mi(ctrE,i) = Data.Analysis{1,i}.MI_i;
%         gdata.CExc.Fr(ctrE,i) = Data.firing_rate(1,i);
%     end
%     gdata.CExc.param(ctrE,1:7) = model_params.param;
%     gdata.CExc.events(ctrE,1:4) = model_params.r_spikevents(2,1:4);
%     gdata.CExc.spikecount(ctrE) =  model_params.r_spikecount;
%     gdata.CExc.eta(ctrE,:) =  model_params.eta;
%     gdata.CExc.gamma(ctrE,:) =  model_params.gamma;
%     gdata.CExc.spikeshape(ctrE,:) =  model_params.spikeshape;
%     gdata.CExc.G2(ctrE) = model_params.G2;
%     gdata.CExc.G5(ctrE) = model_params.G5;
    gdata.Cexc.G10(ctrE) = model_params.G10;
%     gdata.CExc.Mc(ctrE) = model_params.Mc;
    clearvars -except ctrI ctrE CInhfiles CExcfiles gdata D1Inhfiles D1Excfiles
end
% gdata.CExc.param(:,6) = gdata.CExc.param(:,6) + gdata.CExc.Mc';

disp('CExc done')
for ctrE2 = 1:size(D1Excfiles,1)
    load(D1Excfiles(ctrE2).name(1:end-4))
    gdata.D1exc.mspikestimes{ctrE2} = model_params.m_spiketimes*20;
    gdata.D1exc.rspikestimes{ctrE2} = model_params.r_spiketimes*20;
%     for i = 1:size(Data.Analysis,2)
%         gdata.D1Exc.Mi(ctrE2,i) = Data.Analysis{1,i}.MI_i;
%         gdata.D1Exc.Fr(ctrE2,i) = Data.firing_rate(1,i);
%     end
%     gdata.D1Exc.param(ctrE2,1:7) = model_params.param;
%     gdata.D1Exc.events(ctrE2,1:4) = model_params.r_spikevents(2,1:4);
%     gdata.D1Exc.spikecount(ctrE2) =  model_params.r_spikecount;
%     gdata.D1Exc.eta(ctrE2,:) =  model_params.eta;
%     gdata.D1Exc.gamma(ctrE2,:) =  model_params.gamma;
%     gdata.D1Exc.spikeshape(ctrE2,:) =  model_params.spikeshape;
%     gdata.D1Exc.G2(ctrE2) = model_params.G2;
%     gdata.D1Exc.G5(ctrE2) = model_params.G5;
    gdata.D1exc.G10(ctrE2) = model_params.G10;
%     gdata.D1Exc.Mc(ctrE2) = model_params.Mc;
    clearvars -except ctrI ctrE CInhfiles CExcfiles D1Inhfiles D1Excfiles gdata ctrE2 ctrI2
end
% gdata.D1Exc.param(:,6) = gdata.D1Exc.param(:,6) + gdata.D1Exc.Mc';

disp('D1Exc done')
for ctrI2 = 1:size(D1Inhfiles,1)
    load(D1Inhfiles(ctrI2).name(1:end-4))
    gdata.D1inh.mspikestimes{ctrI2} = model_params.m_spiketimes*20;
    gdata.D1inh.rspikestimes{ctrI2} = model_params.r_spiketimes*20;
%     for i = 1:size(Data.Analysis,2)
%         gdata.D1inh.Mi(ctrI2,i) = Data.Analysis{1,i}.MI_i;
%         gdata.D1inh.Fr(ctrI2,i) = Data.firing_rate(1,i);
%     end
%     gdata.D1inh.param(ctrI2,1:7) = model_params.param;
%     gdata.D1inh.events(ctrI2,1:4) = model_params.r_spikevents(2,1:4);
%     gdata.D1inh.spikecount(ctrI2) =  model_params.r_spikecount;
%     gdata.D1inh.eta(ctrI2,:) =  model_params.eta;
%     gdata.D1inh.gamma(ctrI2,:) =  model_params.gamma;
%     gdata.D1inh.spikeshape(ctrI2,:) =  model_params.spikeshape;
%     gdata.D1inh.G2(ctrI2) = model_params.G2;
%     gdata.D1inh.G5(ctrI2) = model_params.G5;
    gdata.D1inh.G10(ctrI2) = model_params.G10;
%     gdata.D1inh.Mc(ctrI2) = model_params.Mc;
    clearvars -except ctrI CInhfiles CExcfiles D1Inhfiles D1Excfiles gdata ctrI2 ctrE2
end
% gdata.D1inh.param(:,6) = gdata.D1inh.param(:,6) + gdata.D1inh.Mc';
disp('D1Inh done')

% [~,gdata.pInhvsExc.param(1)] = ttest2(gdata.Cinh.param(:,1),gdata.CExc.param(:,1));
% [~,gdata.pInhvsExc.param(2)] = ttest2(gdata.Cinh.param(:,2),gdata.CExc.param(:,2));
% [~,gdata.pInhvsExc.param(3)] = ttest2(gdata.Cinh.param(:,3),gdata.CExc.param(:,3));
% [~,gdata.pInhvsExc.param(4)] = ttest2(gdata.Cinh.param(:,4),gdata.CExc.param(:,4));
% [~,gdata.pInhvsExc.param(5)] = ttest2(gdata.Cinh.param(:,5),gdata.CExc.param(:,5));
% [~,gdata.pInhvsExc.param(6)] = ttest2(gdata.Cinh.param(:,6),gdata.CExc.param(:,6));
% [~,gdata.pInhvsExc.param(7)] = ttest2(gdata.Cinh.param(:,7),gdata.CExc.param(:,7));
% [~,gdata.pInhCvsD1.param(1)] = ttest2(gdata.Cinh.param(:,1),gdata.D1inh.param(:,1));
% [~,gdata.pInhCvsD1.param(2)] = ttest2(gdata.Cinh.param(:,2),gdata.D1inh.param(:,2));
% [~,gdata.pInhCvsD1.param(3)] = ttest2(gdata.Cinh.param(:,3),gdata.D1inh.param(:,3));
% [~,gdata.pInhCvsD1.param(4)] = ttest2(gdata.Cinh.param(:,4),gdata.D1inh.param(:,4));
% [~,gdata.pInhCvsD1.param(5)] = ttest2(gdata.Cinh.param(:,5),gdata.D1inh.param(:,5));
% [~,gdata.pInhCvsD1.param(6)] = ttest2(gdata.Cinh.param(:,6),gdata.D1inh.param(:,6));
% [~,gdata.pInhCvsD1.param(7)] = ttest2(gdata.Cinh.param(:,7),gdata.D1inh.param(:,7));
% [~,gdata.pExcCvsD1.param(1)] = ttest2(gdata.CExc.param(:,1),gdata.D1Exc.param(:,1));
% [~,gdata.pExcCvsD1.param(2)] = ttest2(gdata.CExc.param(:,2),gdata.D1Exc.param(:,2));
% [~,gdata.pExcCvsD1.param(3)] = ttest2(gdata.CExc.param(:,3),gdata.D1Exc.param(:,3));
% [~,gdata.pExcCvsD1.param(4)] = ttest2(gdata.CExc.param(:,4),gdata.D1Exc.param(:,4));
% [~,gdata.pExcCvsD1.param(5)] = ttest2(gdata.CExc.param(:,5),gdata.D1Exc.param(:,5));
% [~,gdata.pExcCvsD1.param(6)] = ttest2(gdata.CExc.param(:,6),gdata.D1Exc.param(:,6));
% [~,gdata.pExcCvsD1.param(7)] = ttest2(gdata.CExc.param(:,7),gdata.D1Exc.param(:,7));

keyboard
