function [TotalArea,onset,onsetindex,onsettime,offsettime,onsetcurrent,offsetcurrent,onsetslopecurrent,offsetslopecurrent] = analyze_VCintegral(traceoutput,timetrace,sweepstartindex,protocol)

% A is de first diff of the Trace
A = diff(traceoutput);
Amean = mean(A);
Astd = std(A);

% X will be used for the X-axis of the plot containing the datapoints
% location
X = 1:length(A);

% Below 0 for the onset 
B = X(A < (Amean - (1.5*Astd)));
onsetindex = B(1);
if strcmp(protocol,'100')
    while onsetindex < 500 || B(2) ~= B(1) + 1
        B(1) = [];
        onsetindex = B(1);
    end
elseif strcmp(protocol,'50')
    while onsetindex < 300 || B(2) ~= B(1) + 1 
        B(1) = [];
        onsetindex = B(1);
    end
elseif strcmp(protocol,'10')
    while onsetindex < 75 || B(2) ~= B(1) + 1
        B(1) = [];
        onsetindex = B(1);
    end
end
onsetindex = onsetindex - 3;
onset = traceoutput(onsetindex);

% Above 0 for the offset
C = X(A > (Amean + (0.50*Astd)));
offsetindex = C(find(diff(C) > 5,1));
while offsetindex < onsetindex
    C(1) =[];
    offsetindex = C(find(diff(C) > 5,1));
end
offset = traceoutput(offsetindex);

% Peak between the onset and offset
peakindex = find(traceoutput(onsetindex:offsetindex) == min(traceoutput(onsetindex:offsetindex)),1)+ onsetindex;
peak = traceoutput(peakindex);

% Sloop duration and current amount/ms
onsettime = max(timetrace(onsetindex:peakindex)-timetrace(onsetindex));
offsettime = max(timetrace(peakindex:offsetindex)-timetrace(peakindex));
onsetcurrent = abs(peak-onset);
offsetcurrent = abs(peak-offset);
onsetslopecurrent = (onsetcurrent/onsettime)/1000; % nA/ms
offsetslopecurrent = (offsetcurrent/offsettime)/1000; %nA/ms

% Onset mirror calculation (needed for integrals).
D = find(min(abs(traceoutput(onsetindex+5:offsetindex)-onset)) == abs(traceoutput(onsetindex+5:offsetindex)-onset),1);
onsetmirrorindex = onsetindex + D + 4;
onsetmirror = traceoutput(onsetmirrorindex);

if onsetmirrorindex == offsetindex
    % Total Area
    Area1 = traceoutput(onsetindex:onsetmirrorindex)-onset; % -onset because we want to start at 0
    Area1opp = abs(trapz(timetrace(onsetindex:onsetmirrorindex),Area1));
    TotalArea = Area1opp;
else
    % Area 1
    Area1 = traceoutput(onsetindex:onsetmirrorindex)-onset; % -onset because we want to start at 0
    Area1opp = abs(trapz(timetrace(onsetindex:onsetmirrorindex),Area1));
    
    % Area 2
    Area2 = traceoutput(onsetmirrorindex:offsetindex)-onsetmirror; % -onset because we want to start at 0
    Area2opp = trapz(timetrace(onsetmirrorindex:offsetindex),Area2);
    
    % Area 3
    Yaxes = offset-onset;
    Xaxes = timetrace(offsetindex)-timetrace(onsetindex);
    Area3opp = 0.5*Yaxes*Xaxes;
    
    % TotalAres nA*ms
    TotalArea = Area1opp + (Area3opp - Area2opp);
end

% Check
onsetindex = onsetindex + sweepstartindex;
