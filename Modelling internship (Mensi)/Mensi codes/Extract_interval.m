function [output] = Extract_interval(input,spiketimes,dt)
%
%   Remove the spikes from the data
%   input = vector of real value (i.e. current or voltage)
%   spiketimes = index of the spikes
%   dt = period of time that are removed following each spikes

if(isempty(spiketimes))
    output = input;
else
    output = input(1:spiketimes(1));
    
    for i=1:length(spiketimes)-1
        if(spiketimes(i) + dt <= spiketimes(i+1))
            output = [output;input(spiketimes(i)+dt:spiketimes(i+1))];
        end
    end

    output = [output;input(spiketimes(end)+dt:end)];
end

% OKOKOKOKKOKOKOKOKOKOKOKKOKOKOKOKKOKOOKOKOKOKOKOOKOKOKOKOKOKOOKOKOK %





