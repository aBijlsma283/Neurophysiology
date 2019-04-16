function [] = fit_IF()
%---------------------------------------------------------------------------%
%
%   Extraction and Classification of Three Cortical Neuron Types Reveals Two
%   Distinct Adaptation Mechanisms
%   Skander Mensi, Richard Naud, Christian Pozzorini, Michael Avermann,
%   Carl C.H. Petersen, and Wulfram Gerstner
%   Journal of NeuroPhysiology, 2011
%
%
%   Fit stochastic IF model with spike-triggered current eta and moving
%   threshold gamma as described by Eq 1-3, using a voltage traces V and
%   the injected current I. The fitting procedure is explained in detail
%   in the original paper.
%
%   I = colored noise current in [nA]
%   V = reference voltage generated by injecting I in a IF_eta model [mV]
%   sampling_freq = sampling frequency in [Hz]
%   dt = timestep in [ms]
%
%   There is no free parameters, however one can easily modify the script
%   by adjusting some variables. For instance one can set size_eta to
%   modify the length of the extracted spike-triggered current or one can
%   adjust the number of basis function used to estimate gamma by adjusting
%   nbr_bin_eta
%
%   Reference and Fitted Model:
%   C(V_dot(t)) = -g_l(V(t)-E_l) + I(t) + sum_(t_hat)(eta(t-t_hat))
%   lambda(t) = lambda_0 exp(V(t) - ((V_0 + sum_(t_hat)gamma(t-t_hat))/Delta_V)
%   If spike:   V <- E_reset
%               t = t + t_refr
%               t <- t_hat
%
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
%
%   Main step of the Fitting procedure:
%   1. Extract Spike Shape, T_refr and E_reset from V
%   2. Estimate optimal eta(t) and passive parameters with linear regression
%   3. Estimate moving threshold gamma(t) by maximizing loglikelihood
%   4. Estimate performance of the fitted model on new data set
%   5. Plot results
%
%   Produced figure:
%   1. Spike Shape with T_refr and E_reset
%   2. eta(t) and passive parameters
%   3. gamma(t)
%   4. predicted voltage traces and data voltage traces
%   5. plot results
%
% NOTE: this script can easily be used on in-vitro recodings. When used with 
%       real data. However, there is come important parameters to check:
%       1. be sure that the spikes are well detected
%       2. Double check the estimated t_refr and E_reset
%       3. choose an appropriate value to set the onset of the spikes
%       4. the parameters delay is used two times with different value, be
%       sure to use an appropiate value for each step.
%
% WARNING:  the fitting procedure is fast but requires high memory during the
%           linear regression (matrix inversion).
%           the longer part of the script is the performance evaluation,
%           because one have to generate 1000 spike trains from the reference
%           model and 1000 spike trains from the fitted model to evalutate
%           Md.
%           If to long, one can reduce the duration of the stimuli for the
%           test set (30 seconds) or the number of repetitions (1000) used to
%           compute the PSTH.
%           
%           On a Mac book pro 2.3 GHz with 4 Go RAM, execution time is: 3 minutes
%
% CONTACT:  For any question or remark on the code contact:
%           skander.mensi@epfl.ch
%           
% Copyright Skander Mensi 2011
%
%---------------------------------------------------------------------------%


%0. Generate data-----------------------------------------------------------%
%Initialization
tic
sampling_freq = 20000;                          %[Hz]
dt = 1e3/sampling_freq;                         %[ms]
t_max = 10000/dt;                               %60 seconds in timestep
time = 0:dt:(t_max-1)*dt;                       %for plot
mu_i =0.6+0.6*sin(0.001*time(1:t_max-(100/dt)));%mean of the input current;
sigma_i = 3.5;                                  %standard deviation of the input current
I = mu_i'+sigma_i*randn(t_max-(100/dt),1);      %Gaussian input current
I = smooth([zeros((50/dt),1);...
   I;zeros((50/dt),1)],50); %filter the white noise input current to produce colored noise
% Reference model parameters
ref_C = 0.3; ref_g_l = 0.03; ref_E_l = -65; ref_E_reset = -75;
ref_t_refr = 2; ref_v0 = -45; ref_DeltaV = 1.5;
ref_param = [ref_C ref_g_l ref_E_l ref_E_reset ref_t_refr ref_v0 ref_DeltaV];
ref_eta1 = 0.1; ref_eta_tau1 = 200; ref_eta2 = 0.2; ref_eta_tau2 = 20;
ref_gamma1 = 5; ref_gamma_tau1 = 150; ref_gamma2 = 10; ref_gamma_tau2 = 45;
exp2_func = @(param,x) param(1)*exp(-x/param(2)) + param(3)*exp(-x/param(4));
ref_eta = exp2_func([ref_eta1 ref_eta_tau1 ref_eta2 ref_eta_tau2],time);
ref_gamma = exp2_func([ref_gamma1 ref_gamma_tau1 ref_gamma2 ref_gamma_tau2],time);
V = IF_eta(I,ref_param,ref_eta(1:5000/dt)',ref_gamma(1:5000/dt)',1,sampling_freq); % generate reference voltage traces
'Data generated in: '
toc
%---------------------------------------------------------------------------%


%1. Extract Spike Shape-----------------------------------------------------%
tic
size_spike_shape = 30/dt;                           % size of the spike shape (30 ms)
spiketimes = Extract_spiketimes(V,sampling_freq);   % Extract spiketimes
nbr_spikes = length(spiketimes);
temp_spikes = spiketimes(:,1) - 0.7/dt;             % set the spiketimes 0.7 ms before the maximum of the AP
spike_shape = nan(size_spike_shape,nbr_spikes);     % i.e. just used for plotting
k=1;
for i=1:nbr_spikes-1
    if(temp_spikes(i)+size_spike_shape <= temp_spikes(i+1))
        spike_shape(:,k) = V(temp_spikes(i):temp_spikes(i)+size_spike_shape-1);
        k = k+1;
    end
end
m_spike_shape = nanmean(spike_shape,2);             % Compute the mean spike shape
time_spike_shape = 0:dt:(size_spike_shape-1)*dt;
E_reset = min(m_spike_shape);                       % E_reset is the min of the spike shape
temp = find(m_spike_shape==E_reset);                % t_refr is the time after maximum of the AP
t_refr = temp(end);                                 % when E_reset is reached.
t_refr = t_refr*dt - 0.7 - dt;                      % to be consistent, remove the 0.7ms added for plot
clear temp size_spike_shape temp_spikes k i
'Spike shape extracted in:'
toc
%---------------------------------------------------------------------------%


%2. Extract model parameters with linear regression-------------------------%
% initialize the basis function
tic
size_eta = 2000;                                    % total size of eta in ms
nbr_bink = 40;                                      % nbr bin in eta
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
spike = spiketimes-(delay/dt);             
diff_v =  [diff(V);0]/dt;                       % compute voltage derivative        
diff_v_int = Extract_interval(diff_v,spike,(delay+t_refr)/dt);  % remove spikes from derivative
M = build_M_matrix(V,spike,delay+t_refr,nbr_bink,bink_size,sampling_freq);  % build the matrix that contains the basis function
X = [V ones(length(V),1) I -1*M']';             % build matrix of the regressors
XX = zeros(size(X,1),length(diff_v_int));            % remove spikes from the matrix
for j=1:size(X,1)
	XX(j,:) = Extract_interval(X(j,:)',spike,(delay+t_refr)/dt);
end
b = regress(diff_v_int,XX');                    % linear regression
C = 1/b(3);
g_l = -b(1)*C;
E_l = (b(2)*C)/g_l;
eta = b(4:end)*C;
eta = temp_x'*eta;
time_eta = 0:dt:(length(eta)-1)*dt;
clear X XX diff_v_int M
'Linear Regression performed in:'
toc
%---------------------------------------------------------------------------%


%3. Extract voltage threshold and gamma with Likelihood method--------------%
tic
delay = 0.05;                           % build everything needed to compute the likelihood of a spike train
spike = spiketimes-(delay/dt);     
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
%---cst threshold-----%
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
time_gamma = 0:dt:(length(gamma)-1)*dt;
clear M X XX X_cst XX_cst
'Moving Threshold extracted in:'
toc
'Fitting done'
%---------------------------------------------------------------------------%


%4. Compute Md*, RMSE, raster, PSTH-----------------------------------------%
tic
param = [C g_l E_l E_reset t_refr v0 DeltaV];
error_param = mean((abs(ref_param-param)./ref_param));
% generate test set
t_max = 10000/dt;                               %30 seconds in timestep
time = 0:dt:(t_max-1)*dt;
mu_i =0.4+0.3*sin(0.001*time(1:t_max-(100/dt)));%mean of the input current;
I = mu_i'+sigma_i*randn(t_max-(100/dt),1);
I = smooth([zeros((50/dt),1);I;zeros((50/dt),1)],50);
delta = 4/dt;
nbr_repet = 1000;
nu_d = IF_eta_nu(I,ref_param,ref_eta(1:length(eta))',ref_gamma(1:length(gamma))'...
    ,nbr_repet,sampling_freq);
md = inprod_gamma(nu_d,nu_d,delta);
nu_m = IF_eta_nu(I,param,eta,gamma,nbr_repet,sampling_freq);
mm = inprod_gamma(nu_m,nu_m,delta);
Md = (2*inprod_gamma(nu_d,nu_m,delta))/(md+mm);
nbr_repet = 20;
vv_d = IF_eta(I,ref_param,ref_eta(1:length(eta))',ref_gamma(1:length(gamma))'...
    ,nbr_repet,sampling_freq);
vv_m = IF_eta(I,param,eta,gamma,nbr_repet,sampling_freq);
spiketrain_d = nan(nbr_repet,length(I)); isi_d = [];
spiketrain_m = nan(nbr_repet,length(I)); isi_m = [];
for i=1:nbr_repet
    spikes_d = Extract_spiketimes(vv_d(:,i),sampling_freq);
    spikes_m = Extract_spiketimes(vv_m(:,i),sampling_freq);
    spiketrain_d(i,spikes_d) = i;
    spiketrain_m(i,spikes_m) = -i;
    isi_d = [isi_d [diff(spikes_d);0]'*dt];
    isi_m = [isi_m [diff(spikes_m);0]'*dt];
end
'Perfomance evaluated in:'
toc
'Md ='
Md
'Error on the parameters ='
error_param
%---------------------------------------------------------------------------%

keyboard
%5 Plot and save results----------------------------------------------------%
figure(1),hold on
plot(time_spike_shape,spike_shape,'k')
plot(time_spike_shape,m_spike_shape,'r','linewidth',2)
plot(t_refr+0.7+dt,E_reset+5,'rv','linewidth',2);
axis('tight'),xlabel('time [ms]','Fontsize',16),ylabel('voltage [mV]','Fontsize',16)
lgd = legend('all spikes','mean spike shape','Location','Best');
lgd.FontSize = 14;

figure(2), hold on
subplot(2,5,1),bar(C,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
hold on,plot(1,ref_C,'+r'),ylabel('C [pF]','Fontsize',16),xlim([0,2])
subplot(2,5,2),bar(g_l,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
hold on,plot(1,ref_g_l,'+r'),ylabel('g_l [nS]','Fontsize',16),xlim([0,2])
subplot(2,5,3),bar(E_l,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
hold on,plot(1,ref_E_l,'+r'),ylabel('E_l [mV]','Fontsize',16),xlim([0,2])
subplot(2,5,4),bar(t_refr,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
hold on,plot(1,ref_t_refr,'+r'),ylabel('t_{refr} [ms]','Fontsize',16),xlim([0,2])
subplot(2,5,5),bar(E_reset,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
hold on,plot(1,ref_E_reset,'+r'),ylabel('E_{reset} [mV]','Fontsize',16),xlim([0,2])
subplot(2,5,6),bar(DeltaV,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
hold on,plot(1,ref_DeltaV,'+r'),ylabel('DeltaV [mV]','Fontsize',16),xlim([0,2])
subplot(2,5,7),bar(v0,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
hold on,plot(1,ref_v0,'+r'),ylabel('v0 [mV]','Fontsize',16),xlim([0,2])
legend('fitted parameters','reference parameters','Location','BestOutside')
subplot(2,5,9),bar(Md,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
ylabel('Md [-]','Fontsize',16),xlim([0,2])
subplot(2,5,10),bar(error_param,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
ylabel('error param [-]','Fontsize',16),xlim([0,2])

figure(3),hold on
xlabel('time [ms]','Fontsize',16), ylabel('Spike-triggered current eta [nA]','Fontsize',16)
plot(time_eta,ref_eta(1:length(eta)),'k','linewidth',2)
plot(time_eta,eta,'b','linewidth',2),grid on
lgd = legend('reference eta','fitted eta','Location','Best');
lgd.FontSize = 14;

figure(4),hold on
xlabel('time [ms]','Fontsize',16), ylabel('Moving threshold gamma [mV]','Fontsize',16)
plot(time_gamma,ref_gamma(1:length(gamma)),'k','linewidth',2)
plot(time_gamma,gamma,'g','linewidth',2),grid on
lgd = legend('reference gamma','fitted gamma','Location','Best');
lgd.FontSize = 14;

figure(5),hold on
subplot(3,1,1),plot(time,I,'k'),ylabel('I [nA]','Fontsize',16)
subplot(3,1,2:3),hold on,plot(time,vv_d(:,1),'k'),plot(time,vv_m(:,1),'r')
ylabel('V [mV]','Fontsize',16),xlabel('time [ms]','Fontsize',16)
lgd = legend('reference model','fitted model','Location','Best');
lgd.FontSize = 14;

figure(6),hold on
subplot(5,1,1),plot(time,I,'k'),ylabel('I [nA]','Fontsize',16)
subplot(5,1,2:3),hold on,plot(time,spiketrain_d,'.k'),plot(time,spiketrain_m,'.r')
subplot(5,1,4:5),hold on,plot(time,nu_d,'k'),plot(time,nu_m,'r')
ylabel('PSTH [hz]','Fontsize',16),xlabel('time [ms]','Fontsize',16)
lgd = legend('reference model','fitted model','Location','Best');
lgd.FontSize = 14;
%---------------------------------------------------------------------------%

keyboard
end





