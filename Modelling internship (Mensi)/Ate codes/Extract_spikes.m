function [t_refr, spike_positions, SpikeEvents, Removedspikes] = Extract_spikes(V)
%
%   Extract the poisition and the total count of thespikes in the voltage trace (V). 
%   The voltage trace (V) should be given in mV, samplig frequency (SF) in
%   Hertz.
%
%   Output: T_refr is a scalar which defines the refraction period in ms. Spike_positions is a vector containing the indices 
%   of the spikes. SpikeEvents is a matix with the amount of singlets, duplets, triplets and quartets and their appearance rate.
%   Removedspikes is a scalar which gives the amount of spikes removed.

threshold = 10;                 % Spikes need to be bigger than the threshold to be counted.
spike = 0;                      
zerocross = zeros(length(V),1);
zerocross(V > threshold) = 1; 
ctrskip = 0;

for ctr=1:length(V)
    if ctr > ctrskip
        if zerocross(ctr) == 1
            if spike == 0 || spike ~= 0 && (V(ctr) - min(V(spike_positions(spike):ctr))) > 30
                spike = spike + 1;
                APend = find(zerocross(ctr:end) == 0, 1, 'first') - 2;
                [APmax,APind] = max(V(ctr:ctr+APend));
                spikevalues(spike) = APmax;
                spike_positions(spike) = ctr + APind - 1;
                ctrskip = ctr + APend;
            end
        else
        end
    else
    end
end

LowAmpSpikes = spikevalues < (mean(spikevalues) - 2*std(spikevalues));
spike_positions(LowAmpSpikes) = [];

UnknownSpike = find(diff(spike_positions./20) < 3.5);
Removedspikes = spike_positions(UnknownSpike);
spike_positions(UnknownSpike) = [];

%% Downsampling
realtrain = zeros(7200000,1);
realtrain(spike_positions) = 1;
Nstep = length(realtrain)/3;
newtrain = nan(1,Nstep);
for n = 1:Nstep
    newtrain(n) = sum(realtrain((n-1)*3+1:3*n));
end

%% Autocorellogram
Nspikes = sum(newtrain);
realcorr = xcorr(newtrain,newtrain)./Nspikes;
realcorr = round(realcorr,4);
middle = find(floor(realcorr) == 1);
isi = find(realcorr((middle+1):((middle+1) + 5000)) ~= 0);
t_refr = isi(1)*0.15;
keyboard
if t_refr > 10 % isi is not the same as t_refr. So when there is a big ISI (most of the Exc cells) we change it to a preset 5 ms.
    t_refr = 5;
end
keyboard
%% Burst calculations
% Calculatinos of singlets, duplets, triplets and quartets.

check = 5;
MStimes = diff(spike_positions./20);
duplets = find(MStimes < check*2);
triplets = sum(diff(duplets) == 1);
quartets = sum(diff(find(diff(duplets) == 1)) == 1);
duplets = length(duplets) - triplets - quartets;
triplets = triplets - quartets;
singlets = length(spike_positions) - 2*duplets - 3*triplets - 4*quartets;

eventcount = singlets + duplets + triplets + quartets;
singletrate = singlets/eventcount;
dupletrate = duplets/eventcount;
tripletrate = triplets/eventcount;
quartetrate = quartets/eventcount;
SpikeEvents = [singlets,duplets,triplets,quartets; singletrate,dupletrate,tripletrate,quartetrate];

