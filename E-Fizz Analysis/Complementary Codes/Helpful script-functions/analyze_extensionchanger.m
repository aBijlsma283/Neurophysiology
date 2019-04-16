function analyze_extensionchanger(extension)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl

% Changes the extension of the files in the present folder. Function is
% used to tag specific files with an extension that can be used to seperate
% files in the group analysis.
% nargin == 1 for addition.
% nargin == 0 for deletion.

close all; 
warning off

filesmat = dir('*.mat');

if nargin == 1
    for ctr = 1:length(filesmat)
        [~,name] = fileparts(filesmat(ctr).name);
        newname = [name '_' extension];
        movefile([name '.mat'],[newname '.mat'])    
    end
elseif nargin == 0
    for ctr = 1:length(filesmat)
        [~,name] = fileparts(filesmat(ctr).name);
        parts = strsplit(name,'_');
        parts{end} = {};
        parts(cellfun('isempty',parts)) = [];
        newname = strjoin(parts,'_');
        movefile([name '.mat'],[newname '.mat'])
    end
else 
    disp('nargin has to be 0 (deletion) or 1 (addition)')
end

filespng = dir('*.png');

if nargin == 1
    for ctr = 1:length(filespng)
        [~,name] = fileparts(filespng(ctr).name);
        newname = [name '_' extension];
        movefile([name '.png'],[newname '.png'])    
    end
elseif nargin == 0
    for ctr = 1:length(filespng)
        [~,name] = fileparts(filespng(ctr).name);
        parts = strsplit(name,'_');
        parts{end} = {};
        parts(cellfun('isempty',parts)) = [];
        newname = strjoin(parts,'_');
        movefile([name '.png'],[newname '.png'])
    end
else 
    disp('nargin has to be 0 (deletion) or 1 (addition)')
end

filesfig = dir('*.fig');

if nargin == 1
    for ctr = 1:length(filesfig)
        [~,name] = fileparts(filesfig(ctr).name);
        newname = [name '_' extension];
        movefile([name '.fig'],[newname '.fig'])    
    end
elseif nargin == 0
    for ctr = 1:length(filesfig)
        [~,name] = fileparts(filesfig(ctr).name);
        parts = strsplit(name,'_');
        parts{end} = {};
        parts(cellfun('isempty',parts)) = [];
        newname = strjoin(parts,'_');
        movefile([name '.fig'],[newname '.fig'])
    end
else 
    disp('nargin has to be 0 (deletion) or 1 (addition)')
end
