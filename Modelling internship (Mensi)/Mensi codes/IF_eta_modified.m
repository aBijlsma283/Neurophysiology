function [vv,Mc,spikecount] = IF_eta_modified(I,param,eta,gamma,sampling_freq,nbr_spikes)

dt = 1e3/sampling_freq;
C = param(1); g_l=param(2); E_l=param(3);
E_reset=param(4); t_refr=round(param(5)/dt);
v0 = param(6); DeltaV=param(7);
t_max = length(I);
vv = nan(t_max,1);
Mc = 0;
spikecount = 0;
if isnan(v0) || v0 < 1 || v0 > 1000
    v0 = 30;
    DeltaV = 5;
end

while spikecount < 0.95*nbr_spikes || spikecount > 1.05*nbr_spikes
    if spikecount < 0.95*nbr_spikes
        Mc = Mc - 1;
    elseif spikecount > 1.05*nbr_spikes
        Mc = Mc + 1;
    end
    spike = zeros(1,t_max);
    t_spike = -t_max;
    v=E_l;
    w = zeros(t_max+length(eta),1);
    t=1;
    while(t <= t_max-length(gamma))
        dv = (-g_l*(v-E_l) + I(t) - w(t))/C;
        v = v + dt*dv;
        vt = v0 + Mc;
        vv(round(t),1) = v;
        p = exp(-exp((v-vt)/DeltaV));
        if(rand()>p && t-t_spike>t_refr)
            xmin=-5;
            xmax=5;
            x=xmin+rand(1,1)*(xmax-xmin);
            vv(t,1) = 50 + x;
            vv(t+1:t+t_refr,1) = E_reset;
            t = t+t_refr;
            v = E_reset;
            spike(t) = 1;
            t_spike = t;
            t=t+1;
            w(t:t+length(eta)-1) = w(t:t+length(eta)-1) + eta;
        else
            t=t+1;
        end
    end
    if Mc > 500 || Mc < -500
        error('Mc loop not working')
    end
    
    spikecount = length(nonzeros(vv(:,1) >= 45));
end

