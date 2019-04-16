function [] = FNmodel(filename)
%   Model parameters:
%   C = membrane capacitance
%   g_l = leak conductance
%   E_l = reversal potential
%   E_reset = voltage reset
%   T_refr = absolute refractory period
%   V_0 = voltage threshold baseline
%   DeltaV = stochasticity index
%   eta = dynamics of the spike-triggered current
%   gamma = dynamics of the voltage threshold

%% Loading of the data
warning off
load ([filename '.mat']);
disp(['Fitting: ' filename])

%% Data extraction + starting parameters

V = cell_response * 10e2;
I = input_current * 10e8;
sampling_freq = 20000;                          %[Hz]
dt = 1e3/sampling_freq;                         %[ms]

%% t_refr & Extract Spike Shape
disp('Calculating t_refr and spike events')
[t_refr, model_params.r_spiketimes, model_params.r_spikevents, model_params.r_removedspikes] = Extract_spikes(V);

nbr_spikes = length(model_params.r_spiketimes);
size_spike_shape = 30/dt;                           % size of the spike shape (30 ms)
temp_spikes = (model_params.r_spiketimes(1,:) - 1/dt)';               % set the spiketimes 0.7 ms before the maximum of the AP
spike_shape = nan(size_spike_shape,nbr_spikes);     % i.e. just used for plotting
for i=1:nbr_spikes-1
    if(temp_spikes(i) < length(V)-size_spike_shape)
        spike_shape(:,i) = V(temp_spikes(i):temp_spikes(i)+size_spike_shape-1);
    end
end
m_spike_shape = nanmean(spike_shape,2);             % Compute the mean spike shape
E_reset = min(m_spike_shape);                       % E_reset is the min of the spike shape

clear temp size_spike_shape temp_spikes

%% Extract model parameters with linear regression
disp('Linear regression')
size_eta = 200;                                     % total size of eta in ms
nbr_bink = 30;                                      % nbr bin in eta
temp_bin = @(x,tempdt,bin0) exp(x/tempdt) + bin0;   % function to compute the log-distribution of the kernel bins
tempdt = 0.2;                                       % sharpness of the temp_bin function
x = 0:1/nbr_bink:1-(1/nbr_bink);                    % x-axis of the log-distribution
bink_size0 = 5/dt;                                  % size of the first bin: 5 ms
bin0 = bink_size0*dt;                               % generate the bin distribution and normalize
y = temp_bin(x,tempdt,bin0);
normalization_factor = sum(y); normalization_factor = normalization_factor/size_eta;
bink_size = round(round(y/normalization_factor)/dt);
x=nan(nbr_bink,1);
for j=1:nbr_bink
    if(j==1)
        x(j) = round(bink_size(j)/2.);
    else
        x(j) = round(sum(bink_size(1:j-1))+(bink_size(j)/2.));
    end
end
temp_x = zeros(nbr_bink,sum(bink_size(1:end)));
for j=1:nbr_bink
    temp_x(j,x(j)-round(bink_size(j)/2.)+1:x(j)+round(bink_size(j)/2.)) = ones(1,bink_size(j));
end

% Linear Regression
delay = 1;                                   % parameters used to set the exact timing of the spikes
spike = model_params.r_spiketimes-(delay/dt);
diff_v =  [diff(V);0]/dt;                       % compute voltage derivative
diff_v_int = Extract_interval(diff_v,spike,(delay+t_refr)/dt);  % remove spikes from derivative
M = build_M_matrix(V,spike,delay+t_refr,nbr_bink,bink_size,sampling_freq);  % build the matrix that contains the basis function
X = [V ones(length(V),1) I -1*M']';             % build matrix of the regressors
XX = zeros(size(X,1),length(diff_v_int));            % remove spikes from the matrix
for j=1:size(X,1)
    XX(j,:) = Extract_interval(X(j,:)',spike,(delay+t_refr)/dt);
end
b = pinv(XX')*diff_v_int;                   % linear regression
C = 1/b(3);
g_l = -b(1)*C;
E_l = (b(2)*C)/g_l;
eta = b(4:end)*C;
eta = temp_x'*eta;
time_eta = 0:dt:(length(eta)-1)*dt;
clear X XX M

delay = 0.05;
spike = model_params.r_spiketimes-(delay/dt);
v_int = Extract_interval(V,spike,(delay+t_refr)/dt);
M = build_M_matrix(V,spike,delay+t_refr,nbr_bink,bink_size,sampling_freq);
X = [V -1*ones(length(V),1) -1*M']';
X_cst = [V -1*ones(length(V),1)]';
X_spike = X(:,spike);
X_spike_cst = X_cst(:,spike);
sum_X_spike = sum(X_spike,2)';
sum_X_spike_cst = sum(X_spike_cst,2)';
XX = zeros(size(X,1),length(v_int));
XX_cst = zeros(size(X_cst,1),length(v_int));
for j=1:size(X,1)
    XX(j,:) = Extract_interval(X(j,:)',spike,(delay+t_refr)/dt);
end
for j=1:size(X_cst,1)
    XX_cst(j,:) = Extract_interval(X_cst(j,:)',spike,(delay+t_refr)/dt);
end

% %---cst threshold-----%
theta_0 = [0 0];
tol_theta = 1e-4;
max_iter = 10000;
[theta_cst] = compute_theta(XX_cst,X_spike_cst,sum_X_spike_cst,theta_0,...  % compute v0 for a constant threshold model,
    tol_theta,max_iter);                                                    % used as initial value for the moving threshold
%---dyn threshold-----%

theta_0 = [theta_cst' zeros(1,nbr_bink)];
tol_theta = 1e-4;
max_iter = 50;

[theta] = compute_theta(XX,X_spike,sum_X_spike,theta_0,tol_theta,max_iter); % compute gamma
DeltaV = 1/(theta(1));
v0 = theta(2)*(1/theta(1));
gamma = theta(3:end)*(1/theta(1));
gamma = temp_x'*gamma;
clear M X XX X_cst XX_cst

disp('Model calculations')
param = [C g_l E_l E_reset t_refr v0 DeltaV];
[model_response,Mc] = IF_eta_modified(I,param,eta,gamma,sampling_freq,nbr_spikes);
cell_response = V;
input_current = I;
keyboard
model_params.param = param;
model_params.gamma = gamma;
model_params.eta = eta;
model_params.spikeshape = m_spike_shape;
[~, model_params.m_spiketimes, model_params.m_spikevents, model_params.m_removedspikes] = Extract_spikes(model_response);
model_params.m_spikecount = length(model_params.m_spiketimes);
model_params.r_spikecount = nbr_spikes;
model_params.r_spiketimes = model_params.r_spiketimes./20;
model_params.m_spiketimes = model_params.m_spiketimes./20;
model_params.Mc = Mc;

[~, ~, model_params.G2] = calccofac_symmetric_ignoredoublespikes(model_params.r_spiketimes,model_params.m_spiketimes,model_params.r_spikecount/360000,model_params.m_spikecount/360000,2);
[~, ~, model_params.G5] = calccofac_symmetric_ignoredoublespikes(model_params.r_spiketimes,model_params.m_spiketimes,model_params.r_spikecount/360000,model_params.m_spikecount/360000,5);
[~, ~, model_params.G10] = calccofac_symmetric_ignoredoublespikes(model_params.r_spiketimes,model_params.m_spiketimes,model_params.r_spikecount/360000,model_params.m_spikecount/360000,10);

clearvars -except cell_response hidden_state input_current model_response settings model_params filename

settings.MIanalysis.windowtype = 'dependsontau';
settings.MIanalysis.windowsize = 20000;             % ms, time window for analysis of tau=50 ms files, other tau scaled accordingly.
settings.MIanalysis.factor_ron_roff = 2;            % how much more often in off state than on state
windowtype = 'dependsontau';
windowsize = 20000;

input_current = input_current * 1000;

dt = 1/settings.sampling_rate;

if ~isfield(settings.MIanalysis, 'windowsize')
    windowsize50 = 20000; % 20 s window for analysis at tau = 50ms
else
    windowsize50 = settings.MIanalysis.windowsize;
end
if  settings.tau <= 50
    indexwindow = round(windowsize50/dt);
else
    wfac = settings.tau/50;
    indexwindow = round(wfac*round(windowsize50/dt));
end

ron = 1./(settings.tau*(1+settings.MIanalysis.factor_ron_roff));
roff = settings.MIanalysis.factor_ron_roff*ron;
input_theory = (input_current - (settings.baseline))./ (settings.amplitude_scaling);

Ntime = length(cell_response);
Nwindow = floor(Ntime/indexwindow);

spiketrain = zeros(size(hidden_state));
spiketrain(model_params.r_spiketimes*20) = 1;

for nw = 1:Nwindow
    window_temp = (nw-1)*indexwindow+1:nw*indexwindow;

    hidden_state_temp       = hidden_state(window_temp);
    input_theory_temp       = input_theory(window_temp);
    spiketrain_temp         = spiketrain(window_temp);

    Data.firing_rate{1, nw} = sum(spiketrain_temp)/(indexwindow*dt/1000);

    Data.Analysis{1, nw} = analyze_exp(ron, roff, hidden_state_temp', input_theory_temp, dt, spiketrain_temp');
    Data.Analysis{1, nw}.FI = Data.Analysis{1,nw}.MI/Data.Analysis{1,nw}.MI_i;
    if abs(Data.Analysis{1, nw}.FI)>2
        disp('FI>2, should not be possible')
    end
    if Data.Analysis{1, nw}.MI_i<0
        disp('MI in input negative')    
    end
end

save([filename '_normal'],'cell_response','hidden_state','input_current','model_response','settings','model_params','Data')

% PMT = ones(14,7);
% PMTnames = {'_C_down','_C_up','_g_l_down','_g_l_up','_E_l_down','_E_l_up','_E_reset_down','_E_reset_up','_t_refr_down','_t_refr_up','_v0_down','_v0_up','_DeltaV_down','_DeltaV_up'};
% n = 1;
% for PMTctr = 1:7
%     PMT(n,PMTctr) = 0.5;
%     PMT(n+1,PMTctr) = 2;
%     n = n + 2;
% end
% 
% for ctr = 1:14
%     disp(ctr)
%     [model_response,model_params.m_spikecount] = IF_eta_second(I,param.*PMT(ctr,:),eta,gamma,sampling_freq,Mc);
%     model_params.m_spiketimes = (Extract_spiketimes(model_response,sampling_freq))./20;
%     [~, ~, model_params.G2] = calccofac_symmetric_ignoredoublespikes(model_params.r_spiketimes,model_params.m_spiketimes,model_params.r_spikecount/360000,model_params.m_spikecount/360000,2);
%     [~, ~, model_params.G5] = calccofac_symmetric_ignoredoublespikes(model_params.r_spiketimes,model_params.m_spiketimes,model_params.r_spikecount/360000,model_params.m_spikecount/360000,5);
%     [~, ~, model_params.G10] = calccofac_symmetric_ignoredoublespikes(model_params.r_spiketimes,model_params.m_spiketimes,model_params.r_spikecount/360000,model_params.m_spikecount/360000,10);
%     model_params.param = param.*PMT(ctr,:);
%     name = strcat(filename,PMTnames(ctr));
%     save(name{1},'cell_response','hidden_state','input_current','model_response','settings','model_params')
% end
%  
% 
