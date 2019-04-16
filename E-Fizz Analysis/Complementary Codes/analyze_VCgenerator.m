function [Vmon,Tvec] = analyze_VCgenerator(hold,frequency,totaltime,stimstart,protocolspeed)

% 9-1-2018, Ate Bijlsma, ate.bijlsma@student.ru.nl
% Quick generation of a VC channel for experiments with a missing VC
% channel. Outputs are the channel (Vmon) and a matching time vector.
warning off
% DO NOT MIND THE INTEGER ERROR IN LINE 33/34!

% hold;         Resting membrane potential in mV.
% frequency;    frequency in hZ.
% totaltime;    length of the complete protocol in s.
% stimstart;    Startime of the stimulation in s.
% protocolspeed;Speed of the sawtooth in s.

Vmon = zeros(1,frequency*totaltime); % prelocation 

% Time vector if needed.
Tvec = 1/frequency:1/frequency:1.1;

% Locating the beginning,peak and end of the sawtooth so it is easier to
% fill the sawtooth inbetween these points.
stimbegin = [stimstart; stimstart+2*protocolspeed; stimstart+4*protocolspeed; stimstart+6*protocolspeed; stimstart+8*protocolspeed] .* frequency;
stimpeak = [stimstart+protocolspeed; stimstart+3*protocolspeed; stimstart+5*protocolspeed; stimstart+7*protocolspeed; stimstart+9*protocolspeed] .* frequency;
stimend = (stimstart+10*protocolspeed) * frequency;

% Filling time before and after the stimulus with the hold.
Vmon(1:stimbegin(1)-1) = hold;
Vmon(stimend:end) = hold;

for ctr = 1:5
    if ctr == 5
    Vmon(stimbegin(ctr):stimpeak(ctr)-1) = linspace(hold,-hold,(protocolspeed*frequency));
    Vmon(stimpeak(ctr):stimend-1) = linspace(-hold,hold,(protocolspeed*frequency));       
    else
    Vmon(stimbegin(ctr):stimpeak(ctr)-1) = linspace(hold,-hold,(protocolspeed*frequency));
    Vmon(stimpeak(ctr):stimbegin(ctr+1)-1) = linspace(-hold,hold,(protocolspeed*frequency));
    end
end
end