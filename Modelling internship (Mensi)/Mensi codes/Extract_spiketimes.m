editfunction [spiketimes] = Extract_spiketimes(voltage, sampling_freq)
%
%   Given a voltage trace 'voltage', detect the spiketimes with the upward
%   zero crossing cirterion and set spiketimes to the maximum value
%   following the zero-crossing in a 1 ms windows
%
%   Output: spiketimes is a vector containing the indices of the spikes in
%   timebin
    
    dt = 1e3/sampling_freq;
    k=1;                                                %number of spikes
    voltage_prime = [0;diff(voltage)];                  %voltage derivative
    limite = 10;                                         %zero-crossing criterion
    Dt = round(1/dt);
    limit_t_refr = floor(0.8/dt); 
    t_refr = limit_t_refr + 1;
    spiketimes_t(1,1) = 1;
    
    for i=1:length(voltage)-Dt          %loop over the voltage traces
        if(voltage(i) >= limite && voltage_prime(i) > 0 && t_refr >= limit_t_refr)
            [temp temp_ind] = max(voltage(i:i+Dt));
            minV = min(voltage(spiketimes_t(k):i+temp_ind));
            if temp - minV > 30
                k = k+1;
                t_refr = 0;
                spiketimes_t(k,1) = i+temp_ind;
            else
                t_refr = t_refr + 1;
            end
        else
            t_refr = t_refr + 1;
        end
    end
    spiketimes_t(:,1) = spiketimes_t(:,1)-1;
    spiketimes(:,1) = spiketimes_t(2:end,1);
    
end





