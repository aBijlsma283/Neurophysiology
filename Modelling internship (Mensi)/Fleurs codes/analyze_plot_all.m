all_together_files = dir('*.mat');
for ctr = 1:size(all_together_files,1)
    try
        analyze_plot(all_together_files(ctr).name(1:end-4))
    catch
    end
end