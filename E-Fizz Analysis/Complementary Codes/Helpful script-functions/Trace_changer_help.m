function Trace_changer_help(filename)

load (filename);

files = whos ('Trace*');
output.data.values = eval(files(2).name);
if round(max(output.data.values(:,2)),2) == 0.0700 || round(max(output.data.values(:,2)),2) == 0.0800 || round(max(output.data.values(:,2)),2) == 0.1000 || round(min(output.data.values(:,2)),2) == -0.0800 || round(min(output.data.values(:,2)),2) == -0.0700 || round(min(output.data.values(:,2)),2) == -0.0900
    for Tc = 1:size(files,1)
        Ts = strsplit(files(Tc).name,'_');
        if str2double(Ts{5}) == 3 || str2double(Ts{5}) == 4
            clear(files(Tc).name)
        end
    end
else
    for Tc = 1:4
        Ts = strsplit(files(Tc).name,'_');
        if str2double(Ts{5}) == 1 || str2double(Ts{5}) == 2
            if size(files,1) == 4
                clear(files(1).name)
                clear(files(2).name)
            elseif size(files,1) == 8
                clear(files(1).name)
                clear(files(2).name)
                clear(files(5).name)
                clear(files(6).name)
            elseif size(files,1) == 12
                clear(files(1).name)
                clear(files(2).name)
                clear(files(5).name)
                clear(files(6).name)
                clear(files(9).name)
                clear(files(10).name)
            end
        end
        if str2double(Ts{5}) == 3
            if size(files,1) == 4
                Trace_1_1_1_1 = eval(files(3).name);
                clear(files(3).name);
                Trace_1_1_1_2 = eval(files(4).name);
                clear(files(4).name);
            elseif size(files,1) == 8
                Trace_1_1_1_1 = eval(files(3).name);
                Trace_2_1_1_1 = eval(files(7).name);
                clear(files(3).name);
                clear(files(7).name);
                Trace_1_1_1_2 = eval(files(4).name);
                Trace_2_1_1_2 = eval(files(8).name);
                clear(files(4).name);
                clear(files(8).name);
            elseif size(files,1) == 12
                Trace_1_1_1_1 = eval(files(3).name);
                Trace_2_1_1_1 = eval(files(7).name);
                Trace_3_1_1_1 = eval(files(11).name);
                clear(files(3).name);
                clear(files(7).name);
                clear(files(11).name);
                Trace_1_1_1_2 = eval(files(4).name);
                Trace_2_1_1_2 = eval(files(8).name);
                Trace_3_1_1_2 = eval(files(12).name);
                clear(files(4).name);
                clear(files(8).name);
                clear(files(12).name);
            end
        end
    end
end
clearvars Ts Tc files filename output

save('temp')
