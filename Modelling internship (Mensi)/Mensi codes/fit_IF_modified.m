function [] = fit_IF_modified()
%%   Reference and Fitted Model:
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

%% Data addition---------------------------------------------------------------------------%

files = dir ('*mat');

for ctr2 = 1:size(files,1)
    
    tic
    
    clearvars -except ctr files spike_array;
    file = files(ctr2).name;
    [~,name] = fileparts(files(ctr2).name);
    parts = strsplit(files(ctr2).name,'_');
    load(file);
    disp(['Fitting: ' file])
    
    %     if isfield(settings,'cell_type') == 1
    %         if strcmpi(settings.cell_type,'excitatory')
    %             newname = [name '_Exc'];
    %             movefile([name '.mat'],[newname '.mat'])
    %         elseif strcmpi(settings.cell_type,'pyramidal')
    %             newname = [name '_Exc'];
    %             movefile([name '.mat'],[newname '.mat'])
    %         elseif strcmpi(settings.cell_type,'fast spiking')
    %             newname = [name '_Inh'];
    %             movefile([name '.mat'],[newname '.mat'])
    %         elseif strcmpi(settings.cell_type,'fast_spiking')
    %             newname = [name '_Inh'];
    %             movefile([name '.mat'],[newname '.mat'])
    %         else
    %             disp(['file ' ctr ' error: ' name '. VACANTION TIME!!!!!!!'])
    %         end
    %     else
    %         if settings.tau == 250
    %             newname = [name '_Exc'];
    %             movefile([name '.mat'],[newname '.mat'])
    %         elseif settings.tau == 200
    %             newname = [name '_Exc'];
    %             movefile([name '.mat'],[newname '.mat'])
    %         elseif settings.tau == 50
    %             newname = [name '_Inh'];
    %             movefile([name '.mat'],[newname '.mat'])
    %         end
    %     end
    
    V = cell_response * 10e2;
    I = input_current * 10e8;
    
    %0. Generate data-----------------------------------------------------------%
    %Initialization
    sampling_freq = 20000;                          %[Hz]
    dt = 1e3/sampling_freq;                         %[ms]
    t_max = 10000/dt;                               %60 seconds in timestep
    time = 0:dt:(t_max-1)*dt;                       %for plot
    sigma_i = 3.5;                                  %standard deviation of the input current
    %---------------------------------------------------------------------------%
    
    %1. Extract Spike Shape-----------------------------------------------------%
    size_spike_shape = 30/dt;                           % size of the spike shape (30 ms)
    spiketimes = Extract_spiketimes(V,sampling_freq);   % Extract spiketimes
    nbr_spikes = length(spiketimes);
    spike_array(ctr2) = nbr_spikes;
    temp_spikes = spiketimes(:,1) - 1/dt;               % set the spiketimes 0.7 ms before the maximum of the AP
    spike_shape = nan(size_spike_shape,nbr_spikes);     % i.e. just used for plotting
    k=1; l=1;
    for i=1:nbr_spikes-1
        if(temp_spikes(i)+size_spike_shape <= temp_spikes(i+1))
            spike_shape(:,k) = V(temp_spikes(i):temp_spikes(i)+size_spike_shape-1);
            k = k+1;
        else
            spike_shape_burst(:,l) = V(temp_spikes(i):temp_spikes(i)+size_spike_shape-1);
            l = l+1;
        end
    end
    m_spike_shape = nanmean(spike_shape,2);             % Compute the mean spike shape
    if exist('spike_shape_burst','dir') == 1
        m_spike_shape_burst = nanmean(spike_shape_burst,2); % Compute the mean spike shape
    end
    time_spike_shape = 0:dt:(size_spike_shape-1)*dt;
    E_reset = min(m_spike_shape);                       % E_reset is the min of the spike shape
    temp = find(m_spike_shape==E_reset);                % t_refr is the time after maximum of the AP
    t_refr = temp(end);                                 % when E_reset is reached.
    t_refr = t_refr*dt - 1 - dt;                        % to be consistent, remove the 0.7ms added for plot
    clear temp size_spike_shape temp_spikes
    %---------------------------------------------------------------------------%
    
    %2. Extract model parameters with linear regression-------------------------%
    % initialize the basis function
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
    b = pinv(XX')*diff_v_int;                   % linear regression
    C = 1/b(3);
    g_l = -b(1)*C;
    E_l = (b(2)*C)/g_l;
    eta = b(4:end)*C;
    eta = temp_x'*eta;
    time_eta = 0:dt:(length(eta)-1)*dt;
    clear X XX M
    
    delay = 0.05;             
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
   
    param = [C g_l E_l E_reset t_refr v0 DeltaV];
    [vv_m(:,1),Mc,spikecount] = IF_eta_modified(I,param,eta,gamma,sampling_freq,nbr_spikes);
    disp(ctr2)
    save(['analyzed_' name])
end

toc


