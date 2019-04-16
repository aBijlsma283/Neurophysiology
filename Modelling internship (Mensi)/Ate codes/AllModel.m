modelfiles = dir ('*.mat'); % Retrieves all the .mat files data in the current directory.
for ctr = 1:size(modelfiles,1)
    try
        FNmodel(modelfiles(ctr).name(1:end-4))
    catch
        disp(['error: ' modelfiles(ctr).name(1:end-4)])
        if exist('Errors','dir') == 7
            movefile([modelfiles(ctr).name(1:end-4) '.mat'] ,'Errors')
        else
            mkdir('Errors');
            movefile([modelfiles(ctr).name(1:end-4) '.mat'] ,'Errors')
        end
    end
end
