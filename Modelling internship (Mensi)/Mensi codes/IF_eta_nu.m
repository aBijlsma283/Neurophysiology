function [nu] = IF_eta_nu(I,param,eta,gamma,nbr_repet,sampling_freq)
keyboard
dt = 1e3/sampling_freq;
C = param(1); g_l=param(2); E_l=param(3);
E_reset=param(4); t_refr=param(5)/dt;
v0 = param(6); DeltaV=param(7);
t_max = length(I);
nu = zeros(1,t_max);

for i=1:nbr_repet
    spike = zeros(1,t_max);
    t_spike = -t_max;
    v=E_l;
    w = zeros(t_max+length(eta),1);
    VT = zeros(t_max+length(gamma),1);
    t=1;
    while(t<t_max-t_refr)
        dv = (-g_l*(v-E_l) + I(t) - w(t))/C;
        v = v + dt*dv;
        vt = v0 + VT(t);
        p = exp(-exp((v-vt)/DeltaV));

        if(rand()>p && t-t_spike>t_refr)
            t = t+t_refr;
            v = E_reset;
            spike(t) = 1;
            t_spike = t;
            t=t+1;
            w(t:t+length(eta)-1) = w(t:t+length(eta)-1) + eta;
            VT(t:t+length(gamma)-1) = VT(t:t+length(gamma)-1) + gamma;
        else
            t=t+1;
        end
    end
    nu = nu + spike;
end
nu = nu/nbr_repet;





