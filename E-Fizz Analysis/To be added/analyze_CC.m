function T = analyze_CC (data)
% Analyze_CC is a function in the intracellular e-fizz toolbox.
%
% The function retrieves the data exported from HEKA into MATLAB and
% perform data analysis, assuming that the data is coming from current
% clamp experiments.
%
% Expected organization of the data is:
% Trace_X_Y_Z_W
% X: cells (experiments)
% Y: dataset
% Z: sweeps in each dataset
% W: channels.  In a current clamp experiment there will be two channels
% for each sweep.  The standard configuration in DNP is that the current
% channel is named as 1, and the voltage channel as 2.  If not change this
% assignment in the global variables list at the start of the code.
%
% The only assumptions the code makes is that the minimum height of a spike
% is 30 mV. Otherwise temporal course of the currrent injection, its
% amplitude, response in the voltage channel is calculated on every sweep
% ensuring that if there are experiments that finish without completing the
% entire stimulation protocol the data will still be useful.
%
% Sample entry:
% analyze_CC('HJJP_160523_EtOH_2')

% This code uses a function from signal processing toolbox.

% Tansu Celikel -- celikel@neurophysiology.nl
% V1.   16.06.2016
%
% V1.1  17.06.2016
% What's new in this version:
% 1. Catch current spikes functionality has been added.  This enables
% analysis of data files that are not "successful" in the sense that there
% is not good current clamping.  It comes with a cost that the program can
% no longer accept variable levels of current injection. Comment this
% section 'Catch current spikes' to recover the mentioned functionality.
% 2. It now deals with the variation in current injection pulse length
% even in single dataset.
%
% V1.2  18.06.2016
% What's new in this version:
% 1. To account with various different current clamp protocols, including
% ramp, we now apply a piecewise cubic Hermite interpolating polynomial
% to the current data to extract its envelope.
% 2. Spike count vs Current step size plot now includes a quadratic fit,
% rather than the linear in the previous release.
% 3. Spike threshold is calculated as the peak of second derivative (of the
% Vm).
% 4. Min spike latency and Mean spike latency plots are replaced with first
%
% V2  11.04.2017 (MAJOR RELEASE)
% What's new in this version:
% 1. The code now deals with dual amplifier recordings.  The user can run
% the same code by changing the current and voltage channels in the global
% variables section of the code.
% 2. Various bug fixes and code improvements
% 3. The code now calculates:
% out.exp.AP.latency.window(lp) = out.exp.AP.latency.last(lp)-out.exp.AP.latency.first(lp);
% out.exp.AP.latency.compression(lp) = out.exp.AP.latency.window(lp)/((gI.duration(2)-gI.duration(1))*1000); % fraction of the stimulus window when there is spikes
% spike half width
% absolute and instantenous firing rate
% 4. Segmentation of spikes into different phases to calculate detailed
% spike shape statistics
% 5. Trial data is no longer averaged. -- In a future release, I will
% implement noise measurements, in ancitipation of the frozen nose data
% acquisition.
% 6. Number of figures is increased to two.  Total number of figurines
% plotted is increased to XX.
% 7. Data export into CSV has been inactivated.  Use analyze_CC_group in
% case you need to export data.
%
% V2.1 06.07.2017
% What's new in this version
% 1. Various bug fixes: Earlier versions had assumptions about the data
% organization.  They are eliminated ensuring better error handling.
% 2. AHP code has been modified to handle inhibitory neurons' distinct AHP
% 3. A new figure is added to display the spiking pattern of the neuron as
% a raster and PSTH.
% 4. Current injection detection has been modified to handle those data
% protocols that have the test current injection in the same sweep where
% there is the varying amplitudes of experimental current injections.
%
% V2.1 02.01.2018
% What's new in this version
% 1. AHP detection method is changed.

close all;

load ([data '.mat']);
fl = whos (['Trace*']);

out.file.name = data;
disp (['Current clamp analysis.  Data from: ' data]);

%% global variables
gI.onset = [];
gI.offset = [];
gI.duration = [];
Ichannel = 1;  % 3 if you used the second amplifier in the set-up #1
Vchannel = 2;  % 4 if you used the second amplifier in the set-up #1
minspikeheight = 0.03; % in Volt

% for the AP shape segmentation we will need to have ~20 points
% before the AP peak and ~30 afterwards.  If there are any spikes
% that do not have the necessary data, then remove them from the
% further analysis
tfgs = 50; % window to isolate AP rise
tfgsoff = tfgs+50; % AP decay is slower, so extend the window
%% learn about the experimental conditions, organize the data
ctr = 1;
for lp = 1:size(fl,1)
    
    x = strsplit(fl(lp).name,'_');
    
    for lp2 = 1:size(x,2)
        % info
        out.exp.tag {lp} = x(1);
        out.exp.cell (lp) = str2num(x{2});
        out.exp.dataset (lp) = str2num(x{3});
        out.exp.sweep (lp) = str2num(x{4});
        out.exp.channel (lp) = str2num(x{5});
    end
    
    % group the relevant data from the workspace
    out.data.tracename {lp} = fl(lp).name;
    out.data.values {lp} = eval(fl(lp).name);
    
    out.data.sweepduration (lp) =  out.data.values {lp}(end,1)-out.data.values{lp}(1,1); % in seconds - assumes the first column is time
    out.data.samplingrate (lp) =  round(size(out.data.values {lp},1)/out.data.sweepduration (lp)); % in Hz
    
    % Gather current injection traces into a single matrix
    if out.exp.channel(lp) == Ichannel % current injection
        gI.matrix.raw(ctr,:) = out.data.values{lp}(:,2);
        gI.matrix.dotfilter(ctr,:) = out.data.values{lp}(:,2); % this is to suppress random noise in the I channel -- Set-up #3 shows the spikes in the V channel
        
        if ctr>1
            gI.matrix.dotfilter(ctr,:) = gI.matrix.dotfilter(ctr,:).* gI.matrix.dotfilter(ctr-1,:);
        end
        ctr = ctr+1;
        
    end
end
clear ctr
%% Calculate the time window of current injection

pcd = smooth(mean(gI.matrix.dotfilter),mean(out.data.samplingrate)/100); % smooth with 10 ms running window
[x1,x2,x3,x4] = findpeaks(pcd,'MinPeakHeight',mean(pcd));
%temp_inds = round([x2(1) x2(end)]);
temp_inds = round([x2(1)-100 x2(max(find(diff(x2)<1000))+1)+100]);


%% Calculate response statistics
for lp = 1:length(out.exp.channel)
    disp (num2str(lp))
    if out.exp.channel(lp) == Ichannel % current measurements
        
        out.exp.current.onset.ind(lp) = min(temp_inds); % onset index
        temptime = out.data.values{lp}(:,1);
        out.exp.current.onset.ms(lp) =  temptime (out.exp.current.onset.ind(lp)); % in sec
        gI.onset = out.exp.current.onset.ind(lp);
        
        out.exp.current.offset.ind(lp) = max(temp_inds) ; % offset index
        out.exp.current.offset.ms(lp) = temptime (out.exp.current.offset.ind(lp)); % in sec
        gI.offset = out.exp.current.offset.ind(lp);
        gI.duration = [out.exp.current.onset.ms(lp) out.exp.current.offset.ms(lp)];
        %gI.duration = round(gI.duration*10)./10; % interpolate the number to the closest
        
        out.exp.current.duration(lp) = gI.duration(2)-gI.duration(1); % in sec
        out.exp.current.amplitude(lp) = mean(out.data.values{lp}(out.exp.current.onset.ind(lp):out.exp.current.offset.ind(lp),2))*10^13/10; % in pA
        % note that because the first channel is the current channel, this
        % organization logic will always result out.exp.current variables
        % to have 1 less entry than the voltage channels.  This is not a
        % problem as below we reorganize the data into a matrix, but if you
        % ever call the variables independently remember the missing last
        % index in the current channel, in respect to the voltage.
        %keyboard
        
    elseif out.exp.channel(lp) == Vchannel % voltage measurement
        
        tempsweep = out.data.values{lp};
        
        %tempV = tempsweep(gI.onset-50:gI.offset+50,:); % only the current injection evoked activity
        tempV = tempsweep(gI.onset:gI.offset,:); % only the current injection evoked activity
        
        MPD = out.data.samplingrate (lp)/500; % suppress any peaks within 2 ms of a local maxima
        
        % find spikes
        [pks,locs,w,p] = findpeaks(tempV(:,2), 'MinPeakProminence', minspikeheight,'MinPeakDistance',MPD,'WidthReference','halfheight');
        % to visually confirm the performance of the spike detection,
        % uncomment the following code:
        % findpeaks(tempV(:,2), 'MinPeakProminence', minspikeheight,'WidthReference','halfheight');
        
        
        % for the AP shape segmentation we will need to have ~20 points
        % before the AP peak and ~30 afterwards.  If there are any spikes
        % that do not have the necessary data, then remove them from the
        % further analysis
        if isempty(pks) == 0
            if locs(1)<tfgs+1
                pks = pks(2:end);
                locs = locs(2:end);
                w = w(2:end);
                p = p(2:end);
            end
            
            if isempty(locs) == 0
                if locs(end)+tfgsoff > size(tempV,1)
                    pks = pks(1:end-1);
                    locs = locs(1:end-1);
                    w = w(1:end-1);
                    p = p(1:end-1);
                end
            end
        end
        % store spike timing information
        out.exp.AP.latency.all{lp}=locs;
        % isolate stimulus evoked spikes
        
        if isempty(pks) == 1
            % in case there is no spike in a given trial.  There used to be
            % exception handling in this context but since the indexing
            % method have been adapted to handle changing number of I and V
            % channels, there is no need for it anymore.  In a future
            % release here we will calculate the subthreshold potential
            % statistics.
        else
            disp(['Now processing spike count registry for ' num2str(lp)])
            
            out.exp.AP.count(lp) = length(pks); % Note that due to the current channel there are additional zeros in this array.  This is corrected later on.
            
            % spike latency
            tempsps = (tempV(locs,1)-gI.duration(1))*1000;
            
            out.exp.AP.latency.first(lp) = min(tempsps);  % in ms
            out.exp.AP.latency.last(lp) = max(tempsps); % in ms
            out.exp.AP.latency.window(lp) = out.exp.AP.latency.last(lp)-out.exp.AP.latency.first(lp); % in ms
            out.exp.AP.latency.compression(lp) = out.exp.AP.latency.window(lp)/((gI.duration(2)-gI.duration(1))*1000); % fraction of the stimulus window when there is spikes
            
            % spike interval
            out.exp.AP.ISI.all{lp} = tempsps; % in ms
            tempsps = diff(tempsps);
            %out.exp.AP.latency.meanISI(lp) = mean(tempsps); % in ms
            %out.exp.AP.latency.minISI(lp) = min(tempsps); % in ms
            %out.exp.AP.latency.maxISI(lp) = max(tempsps); % in ms
            %out.exp.AP.latency.medianISI(lp) = median(tempsps); % in ms
            %out.exp.AP.latency.ISIadapt(lp) = ((tempsps(end)-tempsps(1))/tempsps(1))*100; % %
            
            if length(locs) >1 % if there is only 1 spike there is no ISI
                out.exp.AP.ISI.first(lp) = tempsps(1); % in ms
                out.exp.AP.ISI.last(lp) = tempsps(end); % in ms
                out.exp.AP.ISI.mean(lp) = mean(tempsps); % in ms
                out.exp.AP.ISI.min(lp) = min(tempsps); % in ms
                out.exp.AP.ISI.max(lp) = max(tempsps); % in ms
                out.exp.AP.ISI.median(lp) = median(tempsps); % in ms
                out.exp.AP.ISI.adapt(lp) = ((tempsps(end)-tempsps(1))/tempsps(1))*100; % %
            end
            
            % firing rate
            out.exp.AP.FR.abs (lp) = length(tempsps)*(1/((gI.duration(2)-gI.duration(1)))); % in Hz
            out.exp.AP.FR.instantaneous (lp) = 1000/mean(tempsps); % in Hz
            
            % spike half width -- later in the code we also calculate the
            % full width.
            tempsamplingr = (gI.offset-gI.onset)/(max(gI.duration)-min(gI.duration));
            out.exp.AP.halfwidth{lp} = (w./tempsamplingr)*1000; % in ms
            
            % spike amplitude
            out.exp.AP.amplitude{lp} = p; % in raw Vm values - Note in Heleen and Jolien's dataset Vrest is 0.
            
            % out.exp.AP.onset.V {lp} = (pks - p)*1000; % in mV; The voltage spike starts
            
            % spike threshold
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method #1 - faster but there is occasional missing spike
            % [pks2,locs2,w2,p2] = findpeaks(diff(diff(tempV(:,2))),'MinPeakProminence', 0.0005);
            %             if isempty(pks2)==0
            %                 tempspsthrs = tempV(locs2,2)*1000; % in mV
            %
            %                 out.exp.AP.thr.all{lp} = tempspsthrs; % in mV
            %                 out.exp.AP.thr.min(lp) = min(tempspsthrs); % in mV
            %                 out.exp.AP.thr.max(lp) = max(tempspsthrs); % in mV
            %                 out.exp.AP.thr.first(lp) = min(tempspsthrs); % in mV
            %
            %                 if size(pks,1)>1
            %                     out.exp.AP.thr.second(lp) = max(tempspsthrs); % in mV
            %                 end
            %             end
            %
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method #2 -- slower but the threshold is calculated for each
            % spike individually.
            % check whether indexing will trigger an error
            %             try
            %                 tfg = tempV(min(locs)-25:max(locs)+25,2);
            %                 tfgs = 25; % number of sampling points we go back in time from the peak
            %             catch
            %                 tfg = tempV(min(locs)-10:max(locs)+10,2);
            %                 tfgs = 10;
            %             end
            %
            %             % in a future release, here there will be additional codes on
            %             % AP slope detection along with hAP calculations
            %             for lpAP = 1:size(pks,1)
            %                 tempthrV = diff(diff(tempV(locs(lpAP)-tfgs:locs(lpAP)+tfgs,2)));
            %                 lthr.ind (lpAP) = min(find(tempthrV == max(tempthrV))+tfgs+locs(lpAP));
            %
            %                 if size(tempV,1) > lthr.ind(lpAP)
            %                     tempspsthrs (lpAP) = tempV(lthr.ind(lpAP),2)*1000;
            %                 end
            %
            %             end % AP loop
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method #3 -- Introduced in the Version 2 release -- performs
            % full spike shape segmentation -- This works great in the mice
            % but in rat data the rise phase of AP is often best fit with
            % three different slopes.
            %             spdiff = diff(diff(tempV(:,2)));
            %             spthr.mV = (pks-p)*1000;
            %
            %             for lpAP = 1:size(pks,1)
            %
            %                 % isolate the rise phase of each spike :: DEPOLARIZATION
            %                 APrise = smooth(spdiff(locs(lpAP)-tfgs:locs(lpAP)-1),4);
            %
            %                 % second derivative of the peak is the AP threshold  -- as
            %                 % in Huang et al, 2016, 2017
            %                 xp.AP.thr.onset.ind(lpAP) = locs(lpAP)-min(find(APrise==max(APrise)));
            %                 xp.AP.thr.onset.mV(lpAP) = tempV(xp.AP.thr.onset.ind(lpAP),2)*1000; % in mV
            %
            %                 % AP offset when the membrane crosses the same voltage
            %                 % where the AP was triggered :: REPOLARIZATION
            %                 APfall = tempV(locs(lpAP):locs(lpAP)+tfgsoff,2);
            %
            %                 [a,f]=min(abs((APfall*1000)-xp.AP.thr.onset.mV(lpAP)));
            %
            %                 xp.AP.thr.offset.ind (lpAP) = f+1+locs(lpAP); % f+1 ensures that offset is the first Vm recorded past the threshold
            %                 xp.AP.thr.offset.mV (lpAP) = tempV(xp.AP.thr.offset.ind(lpAP),2)*1000; % in mV
            %                 clear a f
            %
            %
            %                 % AHP:: Afterhyperpolarization statistics
            %                 if lpAP> 1 && lpAP < size(pks,1) % AHP count is 1 less than # of APs as the last is either truncated, somatic current is interrupted
            %                     tempAHP = tempV(xp.AP.thr.offset.ind(lpAP-1):xp.AP.thr.onset.ind(lpAP),2);
            %                     [a,f]=min(tempAHP);
            %
            %                     % fall -- hyperpolarization stats
            %                     xp.AP.AHP.peak.ind(lpAP-1) = xp.AP.thr.offset.ind(lpAP-1)+f; % index of AHP peak
            %                     xp.AP.AHP.peak.mV(lpAP-1) = a;% absolute Vm at the AHP peak
            %                     xp.AP.AHP.amplitude(lpAP-1) = a*1000-xp.AP.thr.offset.mV(lpAP-1); % in mV
            %                     xp.AP.AHP.peak.fallduration(lpAP-1) = (f/tempsamplingr)*1000; % in ms
            %                     xp.AP.AHP.peak.fallslope(lpAP-1) = xp.AP.AHP.amplitude(lpAP-1)/xp.AP.AHP.peak.fallduration(lpAP-1); % in mV/ms
            %
            %                     % width and integral
            %                     tempAHPt = tempAHP-tempAHP(1);
            %                     tempAHPt = tempAHPt(2:end);
            %                     [~,g] = min(abs(tempAHPt));
            %
            %                     xp.AP.AHP.width.tailindex(lpAP-1) = xp.AP.thr.offset.ind(lpAP-1)+g; % index of where AHP ends -- it always starts at 1 due to tempAHP calculation
            %                     xp.AP.AHP.width.duration(lpAP-1) = (g/tempsamplingr)*1000; % in ms
            %                     xp.AP.AHP.width.integral(lpAP-1) = sum(tempAHPt(1:g)); % in Volts
            %
            %                     % rise -- redepolarization stats
            %                     xp.AP.AHP.peak.risingduration(lpAP-1) = ((g-f)/tempsamplingr)*1000; % in ms
            %                 xp.AP.AHP.peak.risingslope(lpAP-1) = xp.AP.AHP.amplitude(lpAP-1)/xp.AP.AHP.peak.risingduration(lpAP-1); % in mV/ms -- AHP.amplitude is symettrical for our purposes as its relative height is constant for both the hyper- and depolarization tails
            %
            %             end
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method #4 -- Perform the same segmentation implemented in
            % Method #3 but use the AP height as a token for the AP
            % threshold search.
            
            ntempV = tempV(:,2);
            spdiff = diff(diff(tempV(:,2)));
            spthr.mV = (pks-p)*1000;
            
            for lpAP = 1:size(pks,1)
                
                %%% dV/dt vs mV
                %keyboard
                
                %                 figure; hold on
                %                 for lpAP = 1:size(pks,1)
                %                     gSpike = (ntempV(locs(lpAP)-tfgs:locs(lpAP)+tfgs*2));
                %                     gSpike_d = diff(gSpike);
                %                     gSpike = gSpike (1:end-1); % so that the length is the same with gSpike_d for plotting
                %                     plot(gSpike, gSpike_d)
                %                 end
                %
                
                %figure; plot(tempV(:,2)); hold on
                %plot(xp.AP.thr.onset.ind,xp.AP.thr.onset.mV/1000,'bs')
                %plot(xp.AP.thr.offset.ind,xp.AP.thr.offset.mV/1000,'bs')
                
                %%%% isolate the rise phase of each spike :: DEPOLARIZATION
                APrise = smooth(spdiff(locs(lpAP)-tfgs:locs(lpAP)-1),4);
                
                % second derivative of the peak is the AP threshold  --
                % as in Huang et al, 2016, 2017
                xp.AP.thr.onset.ind(lpAP) = locs(lpAP)-min(find(APrise==max(APrise)));
                % xp.AP.thr.onset.ind(lpAP) = locs(lpAP)-(tfgs-min(find(APrise==max(APrise))));
                
                xp.AP.thr.onset.mV(lpAP) = tempV(xp.AP.thr.onset.ind(lpAP),2)*1000; % in mV
                
                % AP offset when the membrane crosses the same voltage
                % where the AP was triggered :: REPOLARIZATION
                
                APfall = tempV(locs(lpAP):locs(lpAP)+tfgsoff,2);
                
                % the following is the faster solution, but might fail if
                % there are rebounding afterhyperpolarization peaks
                %[a,f]=min(abs((APfall*1000)-xp.AP.thr.onset.mV(lpAP)));
                %xp.AP.thr.offset.ind (lpAP) = f(1)+1+locs(lpAP); % f+1 ensures that offset is the first Vm recorded past the threshold
                
                
                tempfall = APfall-(xp.AP.thr.onset.mV(lpAP)/1000);
                
                [f,~] = find(tempfall<0); % first zero crossing is the first data below the spike threshold
                if isempty(f) % if Vm does not repolarize
                    f = find(tempfall == min(tempfall));
                end
                xp.AP.thr.offset.ind (lpAP) = f(1)+locs(lpAP); % f+1 ensures that offset is the first Vm recorded past the threshold
                xp.AP.thr.offset.mV (lpAP) = tempV(xp.AP.thr.offset.ind(lpAP),2)*1000; % in mV
                clear a f
                
                % AHP:: Afterhyperpolarization statistics
                %if lpAP>1 && lpAP <size(pks,1) % AHP count is 1 less than # of APs as the last is either truncated, somatic current is interrupted
                
                %tempAHP = tempAHP(xp.AP.thr.offset.ind(lpAP));
                
                tempAHP = tempV(:,2)-tempV(xp.AP.thr.offset.ind(lpAP),2); % zero crossing is where the current AP offset is
                tempAHP2 = tempAHP(xp.AP.thr.offset.ind(lpAP):end);
                
                AHPends = min(find(tempAHP2>0));
                
                if isempty (AHPends)
                    AHPends = length(tempAHP2)-1; % The AHP after the last spike in a train might not reduce the Vm down to normalized 0.  This conditional controls for such an exception.
                end
                
                tempAHP = tempAHP(xp.AP.thr.offset.ind(lpAP):xp.AP.thr.offset.ind(lpAP)+AHPends);
                
                
                % tempAHP = tempV(xp.AP.thr.offset.ind(lpAP-1):xp.AP.thr.onset.ind(lpAP),2);
                [a,f]=min(tempAHP);
                
                
                % fall -- hyperpolarization stats
                xp.AP.AHP.peak.ind(lpAP) = xp.AP.thr.offset.ind(lpAP)+f; % index of AHP peak
                xp.AP.AHP.peak.mV(lpAP) = a;% absolute Vm at the AHP peak
                xp.AP.AHP.amplitude(lpAP) = a*1000; % in mV
                xp.AP.AHP.peak.fallduration(lpAP) = (f/tempsamplingr)*1000; % in ms
                xp.AP.AHP.peak.fallslope(lpAP) = xp.AP.AHP.amplitude(lpAP)/xp.AP.AHP.peak.fallduration(lpAP); % in mV/ms
                
                % width and integral
                xp.AP.AHP.width.tailindex(lpAP) = xp.AP.thr.offset.ind(lpAP)+length(tempAHP); % index of where AHP ends -- it always starts at 1 due to tempAHP calculation
                xp.AP.AHP.width.duration(lpAP) = (length(tempAHP)/tempsamplingr)*1000; % in ms
                xp.AP.AHP.width.integral(lpAP) = sum(abs(tempAHP))*1000; % in mV
                
                % rise -- redepolarization stats
                xp.AP.AHP.peak.risingduration(lpAP) = (f/tempsamplingr)*1000; % in ms
                xp.AP.AHP.peak.risingslope(lpAP) = xp.AP.AHP.amplitude(lpAP)/xp.AP.AHP.peak.risingduration(lpAP); % in mV/ms -- AHP.amplitude is symettrical for our purposes as its relative height is constant for both the hyper- and depolarization tails
                
                
                %                     % fall -- hyperpolarization stats
                %                     xp.AP.AHP.peak.ind(lpAP-1) = xp.AP.thr.offset.ind(lpAP-1)+f; % index of AHP peak
                %                     xp.AP.AHP.peak.mV(lpAP-1) = a;% absolute Vm at the AHP peak
                %                     xp.AP.AHP.amplitude(lpAP-1) = a*1000-xp.AP.thr.offset.mV(lpAP-1); % in mV
                %                     xp.AP.AHP.peak.fallduration(lpAP-1) = (f/tempsamplingr)*1000; % in ms
                %                     xp.AP.AHP.peak.fallslope(lpAP-1) = xp.AP.AHP.amplitude(lpAP-1)/xp.AP.AHP.peak.fallduration(lpAP-1); % in mV/ms
                %
                %                     % width and integral
                %                     tempAHPt = tempAHP(f:end)-xp.AP.thr.onset.mV(lpAP)/1000; % in respect to the AP onset threshold of the same spike
                %                     % tempAHPt = tempAHP(f:end)-tempAHP(1); % re. AP offset
                %                     [~,g] = min(abs(tempAHPt));
                %
                %                     xp.AP.AHP.width.tailindex(lpAP-1) = xp.AP.thr.offset.ind(lpAP-1)+g; % index of where AHP ends -- it always starts at 1 due to tempAHP calculation
                %                     xp.AP.AHP.width.duration(lpAP-1) = ((f+g)/tempsamplingr)*1000; % in ms
                %                     xp.AP.AHP.width.integral(lpAP-1) = sum(tempAHP(1:g)); % in mV
                %
                %                     % Uncomment the following line if you would like to
                %                     % visually confirm AHP identification:
                %                     %plot (tempAHPt); title (['Integral max:' num2str(sum(tempAHPt(1:g)))]);
                %                     %keyboard
                %
                %                     % rise -- redepolarization stats
                %                     xp.AP.AHP.peak.risingduration(lpAP-1) = ((g-f)/tempsamplingr)*1000; % in ms
                %                     xp.AP.AHP.peak.risingslope(lpAP-1) = xp.AP.AHP.amplitude(lpAP-1)/xp.AP.AHP.peak.risingduration(lpAP-1); % in mV/ms -- AHP.amplitude is symettrical for our purposes as its relative height is constant for both the hyper- and depolarization tails
                %
                
                
                %end
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tempspsthrs = xp.AP.thr.onset.mV; % this is just so that we can use
            % the same indexing and plotting comments even if we were to change
            % the spike threshold detection methods at some point.
            
            % Uncomment the following code if you would like to visually
            % confirm the performance of the spike segmentation
            %            figure;
            %             plot(tempV(:,2)); hold on;
            %             plot(xp.AP.thr.onset.ind,xp.AP.thr.onset.mV/1000,'ro');
            %             plot(xp.AP.thr.offset.ind,xp.AP.thr.offset.mV/1000,'bs');
            %             %% plot(xp.AP.AHP.width.tailindex,xp.AP.AHP.throffset.mV(2:end-1)/1000,'kd');
            %             axis tight
            %             legend ({'Voltage';'AP threshold';'AP offset';'AHP offset'})
            %
            
            
            out.exp.AP.thr.all{lp} = tempspsthrs; % in mV
            out.exp.AP.thr.min(lp) = min(tempspsthrs); % in mV
            out.exp.AP.thr.mean(lp) = mean(tempspsthrs); % in mV
            out.exp.AP.thr.median(lp) = median(tempspsthrs); % in mV
            out.exp.AP.thr.max(lp) = max(tempspsthrs); % in mV
            out.exp.AP.thr.first(lp) = tempspsthrs(1); % in mV
            out.exp.AP.thr.last(lp) = tempspsthrs(end); % in mV
            out.exp.AP.thr.adaptRate(lp) = ((out.exp.AP.thr.last(lp)-out.exp.AP.thr.first(lp))/out.exp.AP.thr.first(lp))*100; % in percentile - normalized to the first thr
            
            %out.exp.AP.fullwidth{lp} = (tempV(xp.AP.thr.offset.ind,1)-tempV(xp.AP.thr.onset.ind,1))*1000; % in ms
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AHP cataloging
            %%% Not all AHP stats (or the THR, for that matter) are
            %%% exported. If you need anything extra, see xp.AP has already
            %%% the variable of interest listed;  otherwise contact Tansu.
            
            if length(locs)>1 % AHP calculations require a subsequent spike. If there is only 1 spike, then do not calculate AHP
                % latency to peak as calculated from the AP thr crossing in the
                % repolarization phase of the spike
                temp = xp.AP.AHP.peak.fallduration; % in ms
                out.exp.AP.AHP.peaklatency.all{lp} = temp; % in ms
                out.exp.AP.AHP.peaklatency.min(lp) = min(temp);
                out.exp.AP.AHP.peaklatency.mean(lp) = mean(temp);
                out.exp.AP.AHP.peaklatency.max(lp) = max(temp);
                out.exp.AP.AHP.peaklatency.first(lp) = temp(1);
                out.exp.AP.AHP.peaklatency.last(lp) = temp(end);
                out.exp.AP.AHP.peaklatency.adaptRate(lp) = ((out.exp.AP.AHP.peaklatency.last(lp)-out.exp.AP.AHP.peaklatency.first(lp))/out.exp.AP.AHP.peaklatency.first(lp))*100; % in percentile - normalized to the first thr
                
                % Amplitude of the AHP peak, calculated from AP threshold of
                % the corresponding spike
                temp = xp.AP.AHP.amplitude; % in mV
                out.exp.AP.AHP.peakamp.all{lp} = temp; % in mV
                out.exp.AP.AHP.peakamp.min(lp) = min(temp);
                out.exp.AP.AHP.peakamp.mean(lp) = mean(temp);
                out.exp.AP.AHP.peakamp.max(lp) = max(temp);
                out.exp.AP.AHP.peakamp.first(lp) = temp(1);
                out.exp.AP.AHP.peakamp.last(lp) = temp(end);
                out.exp.AP.AHP.peakamp.adaptRate(lp) = ((out.exp.AP.AHP.peakamp.last(lp)-out.exp.AP.AHP.peakamp.first(lp))/out.exp.AP.AHP.peakamp.first(lp))*100; % in percentile - normalized to the first thr
                clear temp
                
                % Duration of the AHP peak, calculated from AP threshold
                % crossing in the repolarization phase of the spike until the
                % AP threshold of the next spike -- note that there will be
                % ~2 less AHP than the total spikes in any given trial as the
                % last spike's AHP is often truncated by the cessation of the
                % current injection.  Note that there could be only N-1 number
                % of AHP for N spikes, in any case.
                temp = xp.AP.AHP.width.duration; % in ms
                out.exp.AP.AHP.duration.all{lp} = temp; % in ms
                out.exp.AP.AHP.duration.min(lp) = min(temp);
                out.exp.AP.AHP.duration.mean(lp) = mean(temp);
                out.exp.AP.AHP.duration.max(lp) = max(temp);
                out.exp.AP.AHP.duration.first(lp) = temp(1);
                out.exp.AP.AHP.duration.last(lp) = temp(end);
                out.exp.AP.AHP.duration.adaptRate(lp) = ((out.exp.AP.AHP.duration.last(lp)-out.exp.AP.AHP.duration.first(lp))/out.exp.AP.AHP.duration.first(lp))*100; % in percentile - normalized to the first thr
                clear temp
                
                % Integral of the AHP.
                temp = xp.AP.AHP.width.integral; % in mV
                out.exp.AP.AHP.integral.all{lp} = temp; % in mV
                out.exp.AP.AHP.integral.min(lp) = min(temp);
                out.exp.AP.AHP.integral.mean(lp) = mean(temp);
                out.exp.AP.AHP.integral.max(lp) = max(temp);
                out.exp.AP.AHP.integral.first(lp) = temp(1);
                out.exp.AP.AHP.integral.last(lp) = temp(end);
                out.exp.AP.AHP.integral.adaptRate(lp) = ((out.exp.AP.AHP.integral.last(lp)-out.exp.AP.AHP.integral.first(lp))/out.exp.AP.AHP.integral.first(lp))*100; % in percentile - normalized to the first thr
                clear temp
                
                % Slope of AHP hyperpolarization
                temp = xp.AP.AHP.peak.fallslope; % in mV/ms
                out.exp.AP.AHP.fallslope.all{lp} = temp; % in mV/ms
                out.exp.AP.AHP.fallslope.min(lp) = min(temp);
                out.exp.AP.AHP.fallslope.mean(lp) = mean(temp);
                out.exp.AP.AHP.fallslope.max(lp) = max(temp);
                out.exp.AP.AHP.fallslope.first(lp) = temp(1);
                out.exp.AP.AHP.fallslope.last(lp) = temp(end);
                out.exp.AP.AHP.fallslope.adaptRate(lp) = ((out.exp.AP.AHP.fallslope.last(lp)-out.exp.AP.AHP.fallslope.first(lp))/out.exp.AP.AHP.fallslope.first(lp))*100; % in percentile - normalized to the first thr
                clear temp
                
                % Slope of AHP depolarization
                temp = xp.AP.AHP.peak.risingslope; % in mV/ms
                out.exp.AP.AHP.risingslope.all{lp} = temp; % in mV/ms
                out.exp.AP.AHP.risingslope.min(lp) = min(temp);
                out.exp.AP.AHP.risingslope.mean(lp) = mean(temp);
                out.exp.AP.AHP.risingslope.max(lp) = max(temp);
                out.exp.AP.AHP.risingslope.first(lp) = temp(1);
                out.exp.AP.AHP.risingslope.last(lp) = temp(end);
                out.exp.AP.AHP.risingslope.adaptRate(lp) = ((out.exp.AP.AHP.risingslope.last(lp)-out.exp.AP.AHP.risingslope.first(lp))/out.exp.AP.AHP.risingslope.first(lp))*100; % in percentile - normalized to the first thr
                clear temp
            end
            
            % An exception: As the current channel is the first, this logic
            % results in uneven number of entries in the data file.  To
            % ease the data analysis, here equalize the vector lengths of
            % the current to the voltage.  Note the value here is
            % irrelevant as it won't be used in any analysis or plotting.
            if lp == size(fl,1)
                out.exp.current.amplitude(lp) = 0;
                out.exp.current.onset.ind(lp) = 0;
                out.exp.current.onset.ms(lp) = 0;
                out.exp.current.offset.ind(lp) = 0;
                out.exp.current.offset.ms(lp) = 0;
                out.exp.current.duration(lp) = 0;
            end
        end
    else
        out.exp.AP.latency.all{lp}=[];
    end
end
%keyboard
%out.exp.AP.latency.all

lst.cell = unique(out.exp.cell);
lst.dataset = unique(out.exp.dataset);
lst.sweep = unique(out.exp.sweep);


%% catch current spikes
% the following code was introduced to account for a few sweeps in the
% dopamine dataset where the currentclamp was not constant.  The caveat is
% that the program can no longer accept variable levels of current
% injection. Comment this section to recover the functionality.

out.exp.current.onset.ind (find(out.exp.current.onset.ind~= median(nonzeros(out.exp.current.onset.ind)))) = median(nonzeros(out.exp.current.onset.ind));
out.exp.current.offset.ind (find(out.exp.current.offset.ind~= median(nonzeros(out.exp.current.offset.ind)))) = median(nonzeros(out.exp.current.offset.ind));

% add zero entries to the end of the out.exp variables so that referencing
% below will work:

out.exp.current.amplitude(size(out.exp.current.amplitude,2)+1:size(out.exp.dataset,2)) = 0;
out.exp.current.duration(size(out.exp.current.duration,2)+1:size(out.exp.dataset,2)) = 0;
out.exp.current.offset.ind(size(out.exp.current.offset.ind,2)+1:size(out.exp.dataset,2)) = out.exp.current.offset.ind(size(out.exp.current.offset.ind,2));
out.exp.current.onset.ind(size(out.exp.current.onset.ind,2)+1:size(out.exp.dataset,2)) = out.exp.current.onset.ind(size(out.exp.current.onset.ind,2));
out.exp.current.offset.ms(size(out.exp.current.offset.ms,2)+1:size(out.exp.dataset,2)) = 0;
out.exp.current.onset.ms(size(out.exp.current.onset.ms,2)+1:size(out.exp.dataset,2)) = 0;

%% organize the data into a tabulated form and export

for lpc = 1:length(lst.cell) % cycle through each cell
    disp (['Processing Cell Number: ' num2str(lst.cell(lpc)) ' ...'])
    temp_ind_cell = find(out.exp.cell==lst.cell(lpc));
    
    for lpd = 1:length(lst.dataset) % organize into repetitions
        disp (['     ... now repetition (i.e. dataset): ' num2str(lst.dataset(lpd)) ])
        temp_ind_rep = find(out.exp.dataset==lst.dataset(lpd));
        
        % sort the current injection amplitude
        tempamp = out.exp.current.amplitude(temp_ind_rep);
        %tempdur = out.exp.current.duration(temp_ind_rep);
        temptail = min(nonzeros(out.exp.current.offset.ind(temp_ind_rep)))-min(nonzeros(out.exp.current.onset.ind(temp_ind_rep)));
        tempchan = out.exp.channel(temp_ind_rep);
        tempcurrents = sort(tempamp(find(tempchan==Ichannel))); % all current steps in this repetition
        tempdata = out.data.tracename(temp_ind_rep);
        
        for lpsteps = 1:length(tempcurrents)
            %disp([num2str(lpsteps)])
            
            tempind = find (tempamp == tempcurrents(lpsteps)); % note that if the user has more than 1 current injection in the same repetition with the same amplitude, the program will break here
            
            if isempty(tempind)==0
                % find the corresponding voltage trace index
                tempvolind = find(strncmp ([tempdata{tempind}(1:end-1) num2str(Vchannel)], out.data.tracename, length([tempdata{tempind}(1:end-1) num2str(Vchannel)])));
                
                %tempvolind = find(strncmp ([out.data.tracename{tempind}(1:end-1) num2str(Vchannel)], out.data.tracename, length(out.data.tracename{tempind})));
                %disp(['index:' num2str(tempvolind)])
                
                out.sorted.cell{lpc}.rep{lpd}.data.columnar.cell(lpsteps) = lst.cell(lpc); % cell no
                out.sorted.cell{lpc}.rep{lpd}.data.columnar.rep(lpsteps) = lst.dataset(lpd); % repetition no
                
                tfts = out.data.values{tempvolind}(out.exp.current.onset.ind(tempind):out.exp.current.offset.ind(tempind),2);
                out.sorted.cell{lpc}.rep{lpd}.data.voltage(1:temptail,lpsteps) = tfts(1:temptail); % voltage trace
                
                out.sorted.cell{lpc}.rep{lpd}.data.columnar.current.amp(lpsteps) = out.exp.current.amplitude(tempind); % in pA
                out.sorted.cell{lpc}.rep{lpd}.data.columnar.current.onset(lpsteps) = out.exp.current.onset.ind(tempind); % index value
                out.sorted.cell{lpc}.rep{lpd}.data.columnar.current.offset(lpsteps) = out.exp.current.offset.ind(tempind); % index value
                out.sorted.cell{lpc}.rep{lpd}.data.columnar.current.duration(lpsteps) = out.exp.current.duration(tempind); % in sec
                
                %                 if exist('out.exp.AP.count')==0
                %                     disp (['There is no AP in ' data '. No further analysis is performed'])
                %                     save (['analyzed_CC_' data '_NO_SPIKES_FOUND.mat'], 'out');
                %                     return
                %
                %                 else
                
                if length(out.exp.AP.ISI.first) >= tempvolind
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APcount(lpsteps) = sum(out.exp.AP.count(tempvolind)); % total nr of spikes
                    %disp(['AP count:' num2str(sum(out.exp.AP.count(tempvolind))) '_at lpsteps:' num2str(lpsteps)])
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.FR.abs(lpsteps) = out.exp.AP.FR.abs(tempvolind); % in Hz; total # of spikes/duration of current injection
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.FR.instantaneous(lpsteps) = out.exp.AP.FR.instantaneous(tempvolind); % in Hz
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APlat.first(lpsteps) =  out.exp.AP.latency.first(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APlat.last(lpsteps) =  out.exp.AP.latency.last(tempvolind);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APlat.window(lpsteps) =  out.exp.AP.latency.window(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APlat.compression(lpsteps) =  out.exp.AP.latency.compression(tempvolind);
                    %end
                    
                    %if exist('out.exp.AP.ISI.first')
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.ISI.first(lpsteps) =  out.exp.AP.ISI.first(tempvolind);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.ISI.last(lpsteps) =  out.exp.AP.ISI.last(tempvolind);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.ISI.min(lpsteps) =  out.exp.AP.ISI.min(tempvolind);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.ISI.max(lpsteps) =  out.exp.AP.ISI.max(tempvolind);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.ISI.mean(lpsteps) =  out.exp.AP.ISI.mean(tempvolind);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.ISI.median(lpsteps) =  out.exp.AP.ISI.median(tempvolind);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.ISI.adaptationrate(lpsteps) =  out.exp.AP.ISI.adapt(tempvolind);
                    %end
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APthr.first(lpsteps) = out.exp.AP.thr.first(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APthr.last(lpsteps) = out.exp.AP.thr.last(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APthr.min(lpsteps) = out.exp.AP.thr.min(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APthr.max(lpsteps) = out.exp.AP.thr.max(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APthr.mean(lpsteps) = out.exp.AP.thr.mean(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APthr.median(lpsteps) = out.exp.AP.thr.median(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APthr.adaptationrate(lpsteps) = out.exp.AP.thr.adaptRate(tempvolind); % in percentile: (last-first)/first*100
                    
                    %if exist('out.exp.AP.ISI.first')
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peaklatency.min(lpsteps)= out.exp.AP.AHP.peaklatency.min(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peaklatency.mean(lpsteps)= out.exp.AP.AHP.peaklatency.mean(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peaklatency.max(lpsteps)= out.exp.AP.AHP.peaklatency.max(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peaklatency.first(lpsteps)= out.exp.AP.AHP.peaklatency.first(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peaklatency.last(lpsteps)= out.exp.AP.AHP.peaklatency.last(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peaklatency.adaptationrate(lpsteps)= out.exp.AP.AHP.peaklatency.adaptRate(tempvolind); % in percentile: (last-first)/first*100
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peakamp.min(lpsteps)= out.exp.AP.AHP.peakamp.min(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peakamp.mean(lpsteps)= out.exp.AP.AHP.peakamp.mean(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peakamp.max(lpsteps)= out.exp.AP.AHP.peakamp.max(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peakamp.first(lpsteps)= out.exp.AP.AHP.peakamp.first(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peakamp.last(lpsteps)= out.exp.AP.AHP.peakamp.last(tempvolind); % in mV
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.peakamp.adaptationrate(lpsteps)= out.exp.AP.AHP.peakamp.adaptRate(tempvolind); % in percentile: (last-first)/first*100
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.duration.min(lpsteps)= out.exp.AP.AHP.duration.min(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.duration.mean(lpsteps)= out.exp.AP.AHP.duration.mean(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.duration.max(lpsteps)= out.exp.AP.AHP.duration.max(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.duration.first(lpsteps)= out.exp.AP.AHP.duration.first(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.duration.last(lpsteps)= out.exp.AP.AHP.duration.last(tempvolind); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.duration.adaptationrate(lpsteps)= out.exp.AP.AHP.duration.adaptRate(tempvolind); % in percentile: (last-first)/first*100
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.integral.min(lpsteps)= out.exp.AP.AHP.integral.min(tempvolind); % in V
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.integral.mean(lpsteps)= out.exp.AP.AHP.integral.mean(tempvolind); % in V
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.integral.max(lpsteps)= out.exp.AP.AHP.integral.max(tempvolind); % in V
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.integral.first(lpsteps)= out.exp.AP.AHP.integral.first(tempvolind); % in V
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.integral.last(lpsteps)= out.exp.AP.AHP.integral.last(tempvolind); % in V
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.integral.adaptationrate(lpsteps)= out.exp.AP.AHP.integral.adaptRate(tempvolind); % in percentile: (last-first)/first*100
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.fallslope.min(lpsteps)= out.exp.AP.AHP.fallslope.min(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.fallslope.mean(lpsteps)= out.exp.AP.AHP.fallslope.mean(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.fallslope.max(lpsteps)= out.exp.AP.AHP.fallslope.max(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.fallslope.first(lpsteps)= out.exp.AP.AHP.fallslope.first(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.fallslope.last(lpsteps)= out.exp.AP.AHP.fallslope.last(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.fallslope.adaptationrate(lpsteps)= out.exp.AP.AHP.fallslope.adaptRate(tempvolind); % in percentile: (last-first)/first*100
                    
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.risingslope.min(lpsteps)= out.exp.AP.AHP.risingslope.min(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.risingslope.mean(lpsteps)= out.exp.AP.AHP.risingslope.mean(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.risingslope.max(lpsteps)= out.exp.AP.AHP.risingslope.max(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.risingslope.first(lpsteps)= out.exp.AP.AHP.risingslope.first(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.risingslope.last(lpsteps)= out.exp.AP.AHP.risingslope.last(tempvolind); % in mV/ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.AHP.risingslope.adaptationrate(lpsteps)= out.exp.AP.AHP.risingslope.adaptRate(tempvolind); % in percentile: (last-first)/first*100
                    %end
                    
                    tempwidth = out.exp.AP.halfwidth{tempvolind};
                    if isempty(tempwidth) ==1 ; tempwidth=0; end % in case there was no spikes
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APhalfwidth.first(lpsteps) =  tempwidth(1); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APhalfwidth.last(lpsteps) =  tempwidth(end);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APhalfwidth.min(lpsteps) =  min(tempwidth);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APhalfwidth.max(lpsteps) =  max(tempwidth);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APhalfwidth.mean(lpsteps) =  mean(tempwidth);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APhalfwidth.median(lpsteps) =  median(tempwidth);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APhalfwidth.adaptationrate(lpsteps) =  ((tempwidth(end)-tempwidth(1))/tempwidth(1))*100; % percent change
                    
                    tempamp2 = out.exp.AP.amplitude{tempvolind};
                    if isempty(tempamp2) ==1 ; tempamp2=0; end % in case there was no spikes
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APamp.first(lpsteps) =  tempamp2(1); % in ms
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APamp.last(lpsteps) =  tempamp2(end);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APamp.min(lpsteps) =  min(tempamp2);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APamp.max(lpsteps) =  max(tempamp2);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APamp.mean(lpsteps) =  mean(tempamp2);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APamp.median(lpsteps) =  median(tempamp2);
                    out.sorted.cell{lpc}.rep{lpd}.data.columnar.APamp.adaptationrate(lpsteps) =  ((tempamp2(end)-tempamp2(1))/tempamp2(1))*100; % percent change
                end
            end
            
        end
        
    end % dataset (repetition) loop
    
    
    % combine data across repetitions
    for lpg = 1:length(out.sorted.cell)
        %keyboard
        
        d = out.sorted.cell{lpg};
        
        % As we are creating matrix representation of the data, if
        % different current amplitudes were applied across different
        % repetitions (which would be against the definition of
        % repetition, by definition), the following analysis would be
        % irrelevant for those experiments.  Then ask Tansu to export
        % single trial data for you.
        
        gdata.columnar.current.amp = [];
        
        gdata.columnar.APcount = [];
        
        gdata.columnar.FR.abs = [];
        gdata.columnar.FR.ins = [];
        
        gdata.columnar.APlat.first = [];
        gdata.columnar.APlat.last = [];
        gdata.columnar.APlat.window = [];
        gdata.columnar.APlat.compression = [];
        
        gdata.columnar.ISI.min = [];
        gdata.columnar.ISI.max = [];
        gdata.columnar.ISI.mean = [];
        gdata.columnar.ISI.median = [];
        gdata.columnar.ISI.AdapRate = [];
        
        gdata.columnar.APthr.first = [];
        gdata.columnar.APthr.last = [];
        gdata.columnar.APthr.min = [];
        gdata.columnar.APthr.max = [];
        gdata.columnar.APthr.mean = [];
        gdata.columnar.APthr.median = [];
        gdata.columnar.APthr.AdapRate = [];
        
        gdata.columnar.APhalfwidth.first = [];
        gdata.columnar.APhalfwidth.last = []; %
        gdata.columnar.APhalfwidth.min = [];
        gdata.columnar.APhalfwidth.max = [];
        gdata.columnar.APhalfwidth.mean = [];
        gdata.columnar.APhalfwidth.median = [];
        gdata.columnar.APhalfwidth.AdapRate = [];
        
        gdata.columnar.APamp.first = [];
        gdata.columnar.APamp.last = [];
        gdata.columnar.APamp.min = [];
        gdata.columnar.APamp.max = [];
        gdata.columnar.APamp.mean = [];
        gdata.columnar.APamp.median = [];
        gdata.columnar.APamp.AdapRate = [];
        
        gdata.columnar.AHP.peaklatency.min = [];
        gdata.columnar.AHP.peaklatency.mean = [];
        gdata.columnar.AHP.peaklatency.max = [];
        gdata.columnar.AHP.peaklatency.first = [];
        gdata.columnar.AHP.peaklatency.last = [];
        gdata.columnar.AHP.peaklatency.AdapRate = [];
        
        gdata.columnar.AHP.peakamp.min = [];
        gdata.columnar.AHP.peakamp.mean = [];
        gdata.columnar.AHP.peakamp.max = [];
        gdata.columnar.AHP.peakamp.first = [];
        gdata.columnar.AHP.peakamp.last = [];
        gdata.columnar.AHP.peakamp.AdapRate = [];
        
        gdata.columnar.AHP.duration.min = [];
        gdata.columnar.AHP.duration.mean = [];
        gdata.columnar.AHP.duration.max = [];
        gdata.columnar.AHP.duration.first = [];
        gdata.columnar.AHP.duration.last = [];
        gdata.columnar.AHP.duration.AdapRate = [];
        
        gdata.columnar.AHP.integral.min = [];
        gdata.columnar.AHP.integral.mean = [];
        gdata.columnar.AHP.integral.max = [];
        gdata.columnar.AHP.integral.first = [];
        gdata.columnar.AHP.integral.last = [];
        gdata.columnar.AHP.integral.AdapRate = [];
        
        gdata.columnar.AHP.fallslope.min = [];
        gdata.columnar.AHP.fallslope.mean = [];
        gdata.columnar.AHP.fallslope.max = [];
        gdata.columnar.AHP.fallslope.first = [];
        gdata.columnar.AHP.fallslope.last = [];
        gdata.columnar.AHP.fallslope.AdapRate = [];
        
        gdata.columnar.AHP.risingslope.min = [];
        gdata.columnar.AHP.risingslope.mean = [];
        gdata.columnar.AHP.risingslope.max = [];
        gdata.columnar.AHP.risingslope.first = [];
        gdata.columnar.AHP.risingslope.last = [];
        gdata.columnar.AHP.risingslope.AdapRate = [];
        
        for lpr = 1:size(d.rep,2)
            
            gdata.columnar.current.amp = [gdata.columnar.current.amp; d.rep{lpr}.data.columnar.current.amp];
            
            gdata.columnar.APcount = [gdata.columnar.APcount; d.rep{lpr}.data.columnar.APcount];
            
            gdata.columnar.FR.abs = [gdata.columnar.FR.abs; d.rep{lpr}.data.columnar.FR.abs];
            gdata.columnar.FR.ins = [gdata.columnar.FR.ins; d.rep{lpr}.data.columnar.FR.instantaneous];
            
            gdata.columnar.APlat.first = [gdata.columnar.APlat.first; d.rep{lpr}.data.columnar.APlat.first];
            gdata.columnar.APlat.last = [gdata.columnar.APlat.last; d.rep{lpr}.data.columnar.APlat.last];
            gdata.columnar.APlat.window = [gdata.columnar.APlat.window; d.rep{lpr}.data.columnar.APlat.window];
            gdata.columnar.APlat.compression = [gdata.columnar.APlat.compression; d.rep{lpr}.data.columnar.APlat.compression];
            
            %if exist('d.rep{lpr}.data.columnar.ISI.max')
            gdata.columnar.ISI.max = [gdata.columnar.ISI.max; d.rep{lpr}.data.columnar.ISI.max];
            gdata.columnar.ISI.mean = [gdata.columnar.ISI.mean; d.rep{lpr}.data.columnar.ISI.mean];
            gdata.columnar.ISI.median = [gdata.columnar.ISI.median; d.rep{lpr}.data.columnar.ISI.median];
            gdata.columnar.ISI.min = [gdata.columnar.ISI.min; d.rep{lpr}.data.columnar.ISI.min];
            gdata.columnar.ISI.AdapRate = [gdata.columnar.ISI.AdapRate; d.rep{lpr}.data.columnar.ISI.adaptationrate];
            %end
            
            gdata.columnar.APthr.first = [gdata.columnar.APthr.first; d.rep{lpr}.data.columnar.APthr.first];
            gdata.columnar.APthr.last = [gdata.columnar.APthr.last; d.rep{lpr}.data.columnar.APthr.last];
            gdata.columnar.APthr.min = [gdata.columnar.APthr.min; d.rep{lpr}.data.columnar.APthr.min];
            gdata.columnar.APthr.max = [gdata.columnar.APthr.max; d.rep{lpr}.data.columnar.APthr.max];
            gdata.columnar.APthr.mean = [gdata.columnar.APthr.mean; d.rep{lpr}.data.columnar.APthr.mean];
            gdata.columnar.APthr.median = [gdata.columnar.APthr.median; d.rep{lpr}.data.columnar.APthr.median];
            gdata.columnar.APthr.AdapRate = [gdata.columnar.APthr.AdapRate; d.rep{lpr}.data.columnar.APthr.adaptationrate];
            
            gdata.columnar.APhalfwidth.first = [gdata.columnar.APhalfwidth.first; d.rep{lpr}.data.columnar.APhalfwidth.first];
            gdata.columnar.APhalfwidth.last = [gdata.columnar.APhalfwidth.last; d.rep{lpr}.data.columnar.APhalfwidth.last];
            gdata.columnar.APhalfwidth.min = [gdata.columnar.APhalfwidth.min; d.rep{lpr}.data.columnar.APhalfwidth.min];
            gdata.columnar.APhalfwidth.max = [gdata.columnar.APhalfwidth.max; d.rep{lpr}.data.columnar.APhalfwidth.max];
            gdata.columnar.APhalfwidth.mean = [gdata.columnar.APhalfwidth.mean; d.rep{lpr}.data.columnar.APhalfwidth.mean];
            gdata.columnar.APhalfwidth.median = [gdata.columnar.APhalfwidth.median; d.rep{lpr}.data.columnar.APhalfwidth.median];
            gdata.columnar.APhalfwidth.AdapRate = [gdata.columnar.APhalfwidth.AdapRate; d.rep{lpr}.data.columnar.APhalfwidth.adaptationrate];
            
            gdata.columnar.APamp.first = [gdata.columnar.APamp.first; d.rep{lpr}.data.columnar.APamp.first];
            gdata.columnar.APamp.last = [gdata.columnar.APamp.last; d.rep{lpr}.data.columnar.APamp.last];
            gdata.columnar.APamp.min = [gdata.columnar.APamp.min; d.rep{lpr}.data.columnar.APamp.min];
            gdata.columnar.APamp.max = [gdata.columnar.APamp.max; d.rep{lpr}.data.columnar.APamp.max];
            gdata.columnar.APamp.mean = [gdata.columnar.APamp.mean; d.rep{lpr}.data.columnar.APamp.mean];
            gdata.columnar.APamp.median = [gdata.columnar.APamp.median; d.rep{lpr}.data.columnar.APamp.median];
            gdata.columnar.APamp.AdapRate = [gdata.columnar.APamp.AdapRate; d.rep{lpr}.data.columnar.APamp.adaptationrate];
            
            %if exist('d.rep{lpr}.data.columnar.AHP.peaklatency.min')
            gdata.columnar.AHP.peaklatency.min = [gdata.columnar.AHP.peaklatency.min; d.rep{lpr}.data.columnar.AHP.peaklatency.min];
            gdata.columnar.AHP.peaklatency.mean = [gdata.columnar.AHP.peaklatency.mean; d.rep{lpr}.data.columnar.AHP.peaklatency.mean];
            gdata.columnar.AHP.peaklatency.max = [gdata.columnar.AHP.peaklatency.max; d.rep{lpr}.data.columnar.AHP.peaklatency.max];
            gdata.columnar.AHP.peaklatency.first = [gdata.columnar.AHP.peaklatency.first; d.rep{lpr}.data.columnar.AHP.peaklatency.first];
            gdata.columnar.AHP.peaklatency.last = [gdata.columnar.AHP.peaklatency.last; d.rep{lpr}.data.columnar.AHP.peaklatency.last];
            gdata.columnar.AHP.peaklatency.AdapRate = [gdata.columnar.AHP.peaklatency.AdapRate; d.rep{lpr}.data.columnar.AHP.peaklatency.adaptationrate];
            
            gdata.columnar.AHP.peakamp.min = [gdata.columnar.AHP.peakamp.min; d.rep{lpr}.data.columnar.AHP.peakamp.min];
            gdata.columnar.AHP.peakamp.mean = [gdata.columnar.AHP.peakamp.mean; d.rep{lpr}.data.columnar.AHP.peakamp.mean];
            gdata.columnar.AHP.peakamp.max = [gdata.columnar.AHP.peakamp.max; d.rep{lpr}.data.columnar.AHP.peakamp.max];
            gdata.columnar.AHP.peakamp.first = [gdata.columnar.AHP.peakamp.first; d.rep{lpr}.data.columnar.AHP.peakamp.first];
            gdata.columnar.AHP.peakamp.last = [gdata.columnar.AHP.peakamp.last; d.rep{lpr}.data.columnar.AHP.peakamp.last];
            gdata.columnar.AHP.peakamp.AdapRate = [gdata.columnar.AHP.peakamp.AdapRate; d.rep{lpr}.data.columnar.AHP.peakamp.adaptationrate];
            
            gdata.columnar.AHP.duration.min = [gdata.columnar.AHP.duration.min; d.rep{lpr}.data.columnar.AHP.duration.min];
            gdata.columnar.AHP.duration.mean = [gdata.columnar.AHP.duration.mean; d.rep{lpr}.data.columnar.AHP.duration.mean];
            gdata.columnar.AHP.duration.max = [gdata.columnar.AHP.duration.max; d.rep{lpr}.data.columnar.AHP.duration.max];
            gdata.columnar.AHP.duration.first = [gdata.columnar.AHP.duration.first; d.rep{lpr}.data.columnar.AHP.duration.first];
            gdata.columnar.AHP.duration.last = [gdata.columnar.AHP.duration.last; d.rep{lpr}.data.columnar.AHP.duration.last];
            gdata.columnar.AHP.duration.AdapRate = [gdata.columnar.AHP.duration.AdapRate; d.rep{lpr}.data.columnar.AHP.duration.adaptationrate];
            
            gdata.columnar.AHP.integral.min = [gdata.columnar.AHP.integral.min; d.rep{lpr}.data.columnar.AHP.integral.min];
            gdata.columnar.AHP.integral.mean = [gdata.columnar.AHP.integral.mean; d.rep{lpr}.data.columnar.AHP.integral.mean];
            gdata.columnar.AHP.integral.max = [gdata.columnar.AHP.integral.max; d.rep{lpr}.data.columnar.AHP.integral.max];
            gdata.columnar.AHP.integral.first = [gdata.columnar.AHP.integral.first; d.rep{lpr}.data.columnar.AHP.integral.first];
            gdata.columnar.AHP.integral.last = [gdata.columnar.AHP.integral.last; d.rep{lpr}.data.columnar.AHP.integral.last];
            gdata.columnar.AHP.integral.AdapRate = [gdata.columnar.AHP.integral.AdapRate; d.rep{lpr}.data.columnar.AHP.integral.adaptationrate];
            
            gdata.columnar.AHP.fallslope.min = [gdata.columnar.AHP.fallslope.min; d.rep{lpr}.data.columnar.AHP.fallslope.min];
            gdata.columnar.AHP.fallslope.mean = [gdata.columnar.AHP.fallslope.mean; d.rep{lpr}.data.columnar.AHP.fallslope.mean];
            gdata.columnar.AHP.fallslope.max = [gdata.columnar.AHP.fallslope.max; d.rep{lpr}.data.columnar.AHP.fallslope.max];
            gdata.columnar.AHP.fallslope.first = [gdata.columnar.AHP.fallslope.first; d.rep{lpr}.data.columnar.AHP.fallslope.first];
            gdata.columnar.AHP.fallslope.last = [gdata.columnar.AHP.fallslope.last; d.rep{lpr}.data.columnar.AHP.fallslope.last];
            gdata.columnar.AHP.fallslope.AdapRate = [gdata.columnar.AHP.fallslope.AdapRate; d.rep{lpr}.data.columnar.AHP.fallslope.adaptationrate];
            
            gdata.columnar.AHP.risingslope.min = [gdata.columnar.AHP.risingslope.min; d.rep{lpr}.data.columnar.AHP.risingslope.min];
            gdata.columnar.AHP.risingslope.mean = [gdata.columnar.AHP.risingslope.mean; d.rep{lpr}.data.columnar.AHP.risingslope.mean];
            gdata.columnar.AHP.risingslope.max = [gdata.columnar.AHP.risingslope.max; d.rep{lpr}.data.columnar.AHP.risingslope.max];
            gdata.columnar.AHP.risingslope.first = [gdata.columnar.AHP.risingslope.first; d.rep{lpr}.data.columnar.AHP.risingslope.first];
            gdata.columnar.AHP.risingslope.last = [gdata.columnar.AHP.risingslope.last; d.rep{lpr}.data.columnar.AHP.risingslope.last];
            gdata.columnar.AHP.risingslope.AdapRate = [gdata.columnar.AHP.risingslope.AdapRate; d.rep{lpr}.data.columnar.AHP.risingslope.adaptationrate];
            %end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot results
        if size(gdata.columnar.current.amp,1) == 1
            xlbls = round(gdata.columnar.current.amp);
        elseif size(gdata.columnar.current.amp,1) >1
            xlbls = round(mean(gdata.columnar.current.amp));
        end
        
        h = figure; hold on
        
        set(h,'PaperOrientation','landscape','PaperType','A4','PaperOrientation','Landscape');
        set(h,'PaperUnits','normalized','PaperPosition', [0 0 1 1]);
        set(h,'Units', 'Normalized', 'Position', [0 0 1 1]);
        set(h,'DefaulttextFontSize',8)
        
        h1 = subplot(5,5,[1:3 6:8]); % spike traces - only from the first repetition
        plot(d.rep{1}.data.voltage(1:end,:)); axis tight; axis off
        title (['Experiment: ' data])
        
        h2 = subplot(5,5,[4 5]); hold on; % spike count
        boxplot(gdata.columnar.APcount,xlbls,'plotstyle','compact')
        plot(smooth(mean(gdata.columnar.APcount,1),'loess'),'r');
        ylabel ('AP count')
        xlabel ('I (pA)')
        % h2.XTickLabel = xlbls;
        
        h3 = subplot(5,5,[9 10]); hold on; % variance in current injection -- should have 0 difference across repetitions
        plot(smooth(mean(gdata.columnar.current.amp,1)),'r', 'LineWidth', 2)
        boxplot(gdata.columnar.current.amp,xlbls,'plotstyle','compact')
        ylabel ('I (pA)')
        xlabel ('Stimulus steps')
        
        % latency plots
        h4 = subplot(5,5,11); hold on; % Latency to first spike
        boxplot(gdata.columnar.APlat.first,xlbls,'plotstyle','compact')
        %plot(smooth(mean(gdata.columnar.APlat.first,1)),'r')
        %xlabel ('Current amplitude (pA)')
        %h4.XTickLabel = xlbls;
        ylabel ('[ms]')
        title ('AP Latency :: 1st AP','FontSize',8)
        h4.XLim = [0 length(xlbls)+1];
        h4.YLim = [0 max(h4.YLim)*1.2];
        
        % before the V1.2, subplot (5,5,12) displayed minimum spike
        % latency:
        % h5 = subplot(5,5,12); hold on; % Minimum spike latency per stimulus intensity
        % boxplot(gdata.columnar.APlat.min,xlbls,'plotstyle','compact')
        % plot(smooth(mean(gdata.columnar.APlat.min,1)),'r')
        %ylabel ('Min spike latency (ms)')
        %xlabel ('Current amplitude (pA)')
        %h5.XTickLabel = xlbls;
        
        h6 = subplot(5,5,12); hold on; % Last AP's latency per stimulus intensity
        boxplot(gdata.columnar.APlat.last,xlbls,'plotstyle','compact')
        %plot(smooth(mean(gdata.columnar.APlat.max,1)),'r')
        %xlabel ('Current amplitude (pA)')
        %h6.XTickLabel = xlbls;
        ylabel ('[ms]')
        title ('AP Latency :: Last AP','FontSize',8)
        h6.XLim = [0 length(xlbls)+1];
        h6.YLim = [0 max(h6.YLim)*1.2];
        
        h61 = subplot(5,5,13); hold on
        boxplot(gdata.columnar.APlat.window,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('AP Latency :: Window','FontSize',8)
        h61.XLim = [0 length(xlbls)+1];
        h61.YLim = [0 max(h61.YLim)*1.1];
        
        h62 = subplot(5,5,14); hold on
        boxplot(gdata.columnar.APlat.compression,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('AP Latency :: Compression','FontSize',8)
        h62.XLim = [0 length(xlbls)+1];
        h62.YLim = [0 max(h62.YLim)*1.1];
        
        h63 = subplot(5,5,15); hold on
        boxplot(gdata.columnar.FR.ins,xlbls,'plotstyle','compact')
        ylabel ('[Hz]')
        title ('Inst. firing rate','FontSize',8)
        h63.XLim = [0 length(xlbls)+1];
        h63.YLim = [0 max(h63.YLim)*1.1];
        
        
        %%%%%%%%%%%%%%%%%%% ISI plots
        h7 = subplot(5,5,16); hold on; % Minimum interspike interval
        boxplot(gdata.columnar.ISI.min,xlbls,'plotstyle','compact')
        %plot(smooth(mean(gdata.columnar.ISI.min,1)),'r')
        %xlabel ('Current amplitude (pA)')
        %h7.XTickLabel = xlbls;
        ylabel ('[ms]')
        title ('ISI :: Min','FontSize',8)
        h7.XLim = [0 length(xlbls)+1];
        h7.YLim = [0 max(h7.YLim)*1.1];
        
        h5 = subplot(5,5,17); hold on;
        boxplot(gdata.columnar.ISI.max,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('ISI :: Max','FontSize',8)
        h5.XLim = [0 length(xlbls)+1];
        h5.YLim = [0 max(h5.YLim)*1.1];
        
        h8 = subplot(5,5,18); hold on;
        boxplot(gdata.columnar.ISI.mean,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('ISI :: Mean','FontSize',8)
        h8.XLim = [0 length(xlbls)+1];
        h8.YLim = [0 max(h8.YLim)*1.1];
        
        h81 = subplot(5,5,19); hold on;
        boxplot(gdata.columnar.ISI.median,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('ISI :: Median','FontSize',8)
        h81.XLim = [0 length(xlbls)+1];
        h81.YLim = [0 max(h81.YLim)*1.1];
        
        h82 = subplot(5,5,20); hold on;
        boxplot(gdata.columnar.ISI.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('ISI :: Adaptation','FontSize',8)
        h82.XLim = [0 length(xlbls)+1];
        h82.YLim = [min(h82.YLim)*1.1 max(h82.YLim)*1.1];
        
        %%%%%%%%%%%%%%%%%%% AP threshold plots
        h9 = subplot(5,5,21); hold on;
        boxplot(gdata.columnar.APthr.first,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('APthr :: 1st AP','FontSize',8)
        h9.XLim = [0 length(xlbls)+1];
        h9.YLim = [min(h9.YLim)*1.1 max(h9.YLim)*1.1];
        xlabel ('I (pA)')
        
        h10 = subplot(5,5,22); hold on;
        boxplot(gdata.columnar.APthr.last,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('APthr :: Last AP','FontSize',8)
        h10.XLim = [0 length(xlbls)+1];
        h10.YLim = [min(h10.YLim)*1.1 max(h10.YLim)*1.1];
        xlabel ('I (pA)')
        
        h11 = subplot(5,5,23); hold on;
        boxplot(gdata.columnar.APthr.mean,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('APthr :: Mean','FontSize',8)
        h11.XLim = [0 length(xlbls)+1];
        h11.YLim = [min(h11.YLim)*1.1 max(h11.YLim)*1.1];
        xlabel ('I (pA)')
        
        h12 = subplot(5,5,24); hold on;
        boxplot(gdata.columnar.APthr.max,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('APthr :: Max','FontSize',8)
        h12.XLim = [0 length(xlbls)+1];
        h12.YLim = [min(h12.YLim)*1.1 max(h12.YLim)*1.1];
        xlabel ('I (pA)')
        
        h14 = subplot(5,5,25); hold on;
        boxplot(gdata.columnar.APthr.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('APthr :: Adaptation','FontSize',8)
        h14.XLim = [0 length(xlbls)+1];
        h14.YLim = [min(h14.YLim)*1.1 max(h14.YLim)*1.1];
        xlabel ('I (pA)')
        
        
        tic; disp('Writing the Figure 1 into file.')
        print(gcf, '-dpdf',  [data '_Cell_' num2str(lst.cell(lpc)),'_CC_Fig1.pdf']);
        print(gcf, '-dpng',  [data '_Cell_' num2str(lst.cell(lpc)),'_CC_Fig1.png']);
        disp('Done!')
        toc
        
        h2f = figure; hold on
        
        set(h2f,'PaperOrientation','landscape','PaperType','A4','PaperOrientation','Landscape');
        set(h2f,'PaperUnits','normalized','PaperPosition', [0 0 1 1]);
        set(h2f,'Units', 'Normalized', 'Position', [0 0 1 1]);
        set(h2f,'DefaulttextFontSize',8)
        
        
        %%%%%%%%%%%%%%%%%%% AP amplitude plots
        h201 = subplot(6,6,1); hold on;
        boxplot(gdata.columnar.APamp.first,xlbls,'plotstyle','compact')
        ylabel ('[V]')
        title ('APamp :: 1st AP','FontSize',8)
        h201.XLim = [0 length(xlbls)+1];
        h201.YLim = [min(h201.YLim)*1.1 max(h201.YLim)*1.1];
        
        h202 = subplot(6,6,2); hold on;
        boxplot(gdata.columnar.APamp.last,xlbls,'plotstyle','compact')
        ylabel ('[V]')
        title ('APamp :: Last AP','FontSize',8)
        h202.XLim = [0 length(xlbls)+1];
        h202.YLim = [min(h202.YLim)*1.1 max(h202.YLim)*1.1];
        
        h203 = subplot(6,6,3); hold on;
        boxplot(gdata.columnar.APamp.mean,xlbls,'plotstyle','compact')
        ylabel ('[V]')
        title ('APamp :: Mean','FontSize',8)
        h203.XLim = [0 length(xlbls)+1];
        h203.YLim = [min(h203.YLim)*1.1 max(h203.YLim)*1.1];
        
        h2031 = subplot(6,6,4); hold on;
        boxplot(gdata.columnar.APamp.min,xlbls,'plotstyle','compact')
        ylabel ('[V]')
        title ('APamp :: Min','FontSize',8)
        h2031.XLim = [0 length(xlbls)+1];
        h2031.YLim = [min(h2031.YLim)*1.1 max(h2031.YLim)*1.1];
        
        h204 = subplot(6,6,5); hold on;
        boxplot(gdata.columnar.APamp.max,xlbls,'plotstyle','compact')
        ylabel ('[V]')
        title ('APamp :: Max','FontSize',8)
        h204.XLim = [0 length(xlbls)+1];
        h204.YLim = [min(h204.YLim)*1.1 max(h204.YLim)*1.1];
        
        h205 = subplot(6,6,6); hold on;
        boxplot(gdata.columnar.APamp.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('APamp :: Adaptation','FontSize',8)
        h205.XLim = [0 length(xlbls)+1];
        h205.YLim = [min(h205.YLim)*1.1 max(h205.YLim)*1.1];
        
        %%%%%%%%%%%%%%%%%%% AP halfwidth plots
        h301 = subplot(6,6,7); hold on;
        boxplot(gdata.columnar.APhalfwidth.first,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('APhalfwidth :: 1st AP','FontSize',8)
        h301.XLim = [0 length(xlbls)+1];
        h301.YLim = [min(h301.YLim)*1.1 max(h301.YLim)*1.1];
        
        h302 = subplot(6,6,8); hold on;
        boxplot(gdata.columnar.APhalfwidth.last,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('APhalfwidth :: Last AP','FontSize',8)
        h302.XLim = [0 length(xlbls)+1];
        h302.YLim = [0 max(h302.YLim)*1.1];
        
        h303 = subplot(6,6,9); hold on;
        boxplot(gdata.columnar.APhalfwidth.mean,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('APhalfwidth :: Mean','FontSize',8)
        h303.XLim = [0 length(xlbls)+1];
        h303.YLim = [0 max(h303.YLim)*1.1];
        
        h3031 = subplot(6,6,10); hold on;
        boxplot(gdata.columnar.APhalfwidth.min,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('APhalfwidth :: Min','FontSize',8)
        h3031.XLim = [0 length(xlbls)+1];
        h3031.YLim = [0 max(h3031.YLim)*1.1];
        
        h304 = subplot(6,6,11); hold on;
        boxplot(gdata.columnar.APhalfwidth.max,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('APhalfwidth :: Max','FontSize',8)
        h304.XLim = [0 length(xlbls)+1];
        h304.YLim = [0 max(h304.YLim)*1.1];
        
        h305 = subplot(6,6,12); hold on;
        boxplot(gdata.columnar.APhalfwidth.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('APhalfwidth :: Adaptation','FontSize',8)
        h305.XLim = [0 length(xlbls)+1];
        h305.YLim = [min(h305.YLim)*1.1 max(h305.YLim)*1.1];
        
        
        %%%%%%%%%%%%%%%%%%% AHP -- latency to peak
        h401 = subplot(6,6,13); hold on;
        boxplot(gdata.columnar.AHP.peaklatency.first,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('AHP lat2peak :: 1st AP','FontSize',8)
        h401.XLim = [0 length(xlbls)+1];
        h401.YLim = [min(h401.YLim)*1.1 max(h401.YLim)*1.1];
        
        h402 = subplot(6,6,19); hold on;
        boxplot(gdata.columnar.AHP.peaklatency.last,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('AHP lat2peak :: Last AP','FontSize',8)
        h402.XLim = [0 length(xlbls)+1];
        h402.YLim = [min(h402.YLim)*1.1 max(h402.YLim)*1.1];
        
        h502 = subplot(6,6,25); hold on;
        boxplot(gdata.columnar.AHP.peaklatency.mean,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('AHP lat2peak :: Mean','FontSize',8)
        h502.XLim = [0 length(xlbls)+1];
        h502.YLim = [min(h502.YLim)*1.1 max(h502.YLim)*1.1];
        
        h602 = subplot(6,6,31); hold on;
        boxplot(gdata.columnar.AHP.peaklatency.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('AHP lat2peak :: Adaptation','FontSize',8)
        h602.XLim = [0 length(xlbls)+1];
        h602.YLim = [min(h602.YLim)*1.1 max(h602.YLim)*1.1];
        xlabel ('I (pA)')
        
        %%%%%%%%%%%%%%%%%%% AHP -- peak amplitude
        h4011 = subplot(6,6,14); hold on;
        boxplot(gdata.columnar.AHP.peakamp.first,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('AHP peakAmp :: 1st AP','FontSize',8)
        h4011.XLim = [0 length(xlbls)+1];
        h4011.YLim = [min(h4011.YLim)*1.1 max(h4011.YLim)*1.1];
        
        h4021 = subplot(6,6,20); hold on;
        boxplot(gdata.columnar.AHP.peakamp.last,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('AHP peakAmp :: Last AP','FontSize',8)
        h4021.XLim = [0 length(xlbls)+1];
        h4021.YLim = [min(h4021.YLim)*1.1 max(h4021.YLim)*1.1];
        
        h5021 = subplot(6,6,26); hold on;
        boxplot(gdata.columnar.AHP.peakamp.mean,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('AHP peakAmp :: Mean','FontSize',8)
        h5021.XLim = [0 length(xlbls)+1];
        h5021.YLim = [min(h5021.YLim)*1.1 max(h5021.YLim)*1.1];
        
        h6021 = subplot(6,6,32); hold on;
        boxplot(gdata.columnar.AHP.peakamp.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('AHP peakAmp :: Adaptation','FontSize',8)
        h6021.XLim = [0 length(xlbls)+1];
        h6021.YLim = [min(h6021.YLim)*1.1 max(h6021.YLim)*1.1];
        xlabel ('I (pA)')
        
        %%%%%%%%%%%%%%%%%%% AHP -- duration
        h4012 = subplot(6,6,15); hold on;
        boxplot(gdata.columnar.AHP.duration.first,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('AHP duration :: 1st AP','FontSize',8)
        h4012.XLim = [0 length(xlbls)+1];
        h4012.YLim = [min(h4012.YLim)*1.1 max(h4012.YLim)*1.1];
        
        h4022 = subplot(6,6,21); hold on;
        boxplot(gdata.columnar.AHP.duration.last,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('AHP duration :: Last AP','FontSize',8)
        h4022.XLim = [0 length(xlbls)+1];
        h4022.YLim = [min(h4022.YLim)*1.1 max(h4022.YLim)*1.1];
        
        h5022 = subplot(6,6,27); hold on;
        boxplot(gdata.columnar.AHP.duration.mean,xlbls,'plotstyle','compact')
        ylabel ('[ms]')
        title ('AHP duration :: Mean','FontSize',8)
        h5022.XLim = [0 length(xlbls)+1];
        h5022.YLim = [min(h5022.YLim)*1.1 max(h5022.YLim)*1.1];
        
        h6022 = subplot(6,6,33); hold on;
        boxplot(gdata.columnar.AHP.duration.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('AHP duration :: Adaptation','FontSize',8)
        h6022.XLim = [0 length(xlbls)+1];
        h6022.YLim = [min(h6022.YLim)*1.1 max(h6022.YLim)*1.1];
        xlabel ('I (pA)')
        
        %%%%%%%%%%%%%%%%%%% AHP -- integral, area under curve (total
        %%%%%%%%%%%%%%%%%%% voltage)
        h4013 = subplot(6,6,16); hold on;
        boxplot(gdata.columnar.AHP.integral.first,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('AHP integral :: 1st AP','FontSize',8)
        h4013.XLim = [0 length(xlbls)+1];
        h4013.YLim = [min(h4013.YLim)*1.1 max(h4013.YLim)*1.1];
        
        h4023 = subplot(6,6,22); hold on;
        boxplot(gdata.columnar.AHP.integral.last,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('AHP integral :: Last AP','FontSize',8)
        h4023.XLim = [0 length(xlbls)+1];
        h4023.YLim = [min(h4023.YLim)*1.1 max(h4023.YLim)*1.1];
        
        h5023 = subplot(6,6,28); hold on;
        boxplot(gdata.columnar.AHP.integral.mean,xlbls,'plotstyle','compact')
        ylabel ('[mV]')
        title ('AHP integral :: Mean','FontSize',8)
        h5023.XLim = [0 length(xlbls)+1];
        h5023.YLim = [min(h5023.YLim)*1.1 max(h5023.YLim)*1.1];
        
        h6023 = subplot(6,6,34); hold on;
        boxplot(gdata.columnar.AHP.integral.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('AHP integral :: Adaptation','FontSize',8)
        h6023.XLim = [0 length(xlbls)+1];
        h6023.YLim = [min(h6023.YLim)*1.1 max(h6023.YLim)*1.1];
        xlabel ('I (pA)')
        
        %%%%%%%%%%%%%%%%%%% AHP -- fallslope (mV/ms)
        h4014 = subplot(6,6,17); hold on;
        boxplot(gdata.columnar.AHP.fallslope.first,xlbls,'plotstyle','compact')
        ylabel ('[mV/ms]')
        title ('AHP fallslope :: 1st AP','FontSize',8)
        h4014.XLim = [0 length(xlbls)+1];
        h4014.YLim = [min(h4014.YLim)*1.1 max(h4014.YLim)*1.1];
        
        h4024 = subplot(6,6,23); hold on;
        boxplot(gdata.columnar.AHP.fallslope.last,xlbls,'plotstyle','compact')
        ylabel ('[mV/ms]')
        title ('AHP fallslope :: Last AP','FontSize',8)
        h4024.XLim = [0 length(xlbls)+1];
        h4024.YLim = [min(h4024.YLim)*1.1 max(h4024.YLim)*1.1];
        
        h5024 = subplot(6,6,29); hold on;
        boxplot(gdata.columnar.AHP.fallslope.mean,xlbls,'plotstyle','compact')
        ylabel ('[mV/ms]')
        title ('AHP fallslope :: Mean','FontSize',8)
        h5024.XLim = [0 length(xlbls)+1];
        h5024.YLim = [min(h5024.YLim)*1.1 max(h5024.YLim)*1.1];
        
        h6024 = subplot(6,6,35); hold on;
        boxplot(gdata.columnar.AHP.fallslope.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('AHP fallslope :: Adaptation','FontSize',8)
        h6024.XLim = [0 length(xlbls)+1];
        h6024.YLim = [min(h6024.YLim)*1.1 max(h6024.YLim)*1.1];
        xlabel ('I (pA)')
        
        %%%%%%%%%%%%%%%%%%% AHP -- riseslope (mV/ms)
        h4015 = subplot(6,6,18); hold on;
        boxplot(gdata.columnar.AHP.risingslope.first,xlbls,'plotstyle','compact')
        ylabel ('[mV/ms]')
        title ('AHP riseslope :: 1st AP','FontSize',8)
        h4015.XLim = [0 length(xlbls)+1];
        h4015.YLim = [min(h4015.YLim)*1.1 max(h4015.YLim)*1.1];
        
        h4025 = subplot(6,6,24); hold on;
        boxplot(gdata.columnar.AHP.risingslope.last,xlbls,'plotstyle','compact')
        ylabel ('[mV/ms]')
        title ('AHP riseslope :: Last AP','FontSize',8)
        h4025.XLim = [0 length(xlbls)+1];
        h4025.YLim = [min(h4025.YLim)*1.1 max(h4025.YLim)*1.1];
        
        h5025 = subplot(6,6,30); hold on;
        boxplot(gdata.columnar.AHP.risingslope.mean,xlbls,'plotstyle','compact')
        ylabel ('[mV/ms]')
        title ('AHP riseslope :: Mean','FontSize',8)
        h5025.XLim = [0 length(xlbls)+1];
        h5025.YLim = [min(h5025.YLim)*1.1 max(h5025.YLim)*1.1];
        
        h6025 = subplot(6,6,36); hold on;
        boxplot(gdata.columnar.AHP.risingslope.AdapRate,xlbls,'plotstyle','compact')
        ylabel ('[%]')
        title ('AHP riseslope :: Adaptation','FontSize',8)
        h6025.XLim = [0 length(xlbls)+1];
        h6025.YLim = [min(h6025.YLim)*1.1 max(h6025.YLim)*1.1];
        xlabel ('I (pA)')
        
        tic; disp('Writing the Figure 2 into file.')
        print(gcf, '-dpdf',  [data '_Cell_' num2str(lst.cell(lpc)),'_CC_Fig2.pdf']);
        print(gcf, '-dpng',  [data '_Cell_' num2str(lst.cell(lpc)),'_CC_Fig2.png']);
        disp('Done!');
        toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spike plots
        %%% to be included in future release: Phase plots
        
        %keyboard
        h = figure; hold on
        
        set(h,'PaperOrientation','portrait','PaperType','A4','PaperOrientation','Landscape');
        set(h,'PaperUnits','normalized','PaperPosition', [0 0 1 1]);
        set(h,'Units', 'Normalized', 'Position', [0 0 1 1]);
        set(h,'DefaulttextFontSize',8)
        
        %%% create spike matrix - for PSTH
        gspikes = zeros(size(out.exp.AP.latency.all,2),max(out.exp.current.offset.ind)-min(out.exp.current.onset.ind));
        [~,strials] = sort(out.exp.sweep);
        for lpsm = 1:size(out.exp.AP.latency.all,2)
            sind = strials (lpsm);
            if isempty(out.exp.AP.latency.all{sind})==0
                stempind = out.exp.AP.latency.all{sind}';
                gspikes(sind,stempind)=1;
                
                subplot(3,1,[2 3]); hold on
                plot (stempind,lpsm*ones(1,length(stempind)),'k.','MarkerSize',20)
                plot (stempind,lpsm*ones(1,length(stempind)),'ro','MarkerSize',10)
                
            end
        end
        xlabel ('Time')
        ylabel ('Trials (Iamp ranked, and grouped)')
        set(gca,'YLim',[0 size(out.exp.AP.latency.all,2)+1],'Box','on')
        
        subplot(3,1,1)
        plot(sum(gspikes),'k')
        ylabel ('#')
        
        tic; disp('Writing the Figure 3 into file.')
        print(gcf, '-dpdf',  [data '_Cell_' num2str(lst.cell(lpc)),'_CC_Fig3.pdf']);
        print(gcf, '-dpng',  [data '_Cell_' num2str(lst.cell(lpc)),'_CC_Fig3.png']);
        disp('Done!')
        
        toc
        
        %%% in versions prior to V2, the following code exported data in
        %%% CSV format for further data analysis.  The new analyze_CC_group
        %%% function perform the same job and enable combining data across
        %%% different experiments.  Thus the following has been inactivated
        %%% until there is a specific use for it.
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% export data for external analysis
        %         T = table (mean(gdata.columnar.current.amp,1)', mean(gdata.columnar.APcount,1)', ...
        %             mean(gdata.columnar.APlat.first,1)',mean(gdata.columnar.APlat.min,1)', ...
        %             mean(gdata.columnar.APlat.max,1)', mean(gdata.columnar.APlat.mean,1)', ...
        %             mean(gdata.columnar.APthr.first,1)', mean(gdata.columnar.APthr.max,1)', ...
        %             mean(gdata.columnar.APamp.first,1)', mean(gdata.columnar.APamp.min,1)', ...
        %             mean(gdata.columnar.APamp.mean,1)', mean(gdata.columnar.APamp.AdapRate,1)', ...
        %             mean(gdata.columnar.APwidth.first,1)', mean(gdata.columnar.APwidth.min,1)', ...
        %             mean(gdata.columnar.APwidth.mean,1)', mean(gdata.columnar.APwidth.AdapRate,1)', 'VariableNames', ...
        %             {'Current', 'AP_count', 'Latency2firstSpike', 'minAPlatency', ...
        %             'maxAPlatency', 'meanAPlatency', 'APthr_First', 'APthr_Max', 'FirstSpikeAmp', 'minSpikeAmp', ...
        %             'meanSpikeAmp','SpikeAmpAdaptation', 'FirstSpikeWidth', ...
        %             'minSpikeWidth', 'meanSpikeWidth', 'SpikeWidthAdaptation' });
        %
        %         writetable(T,[data '_Cell_' num2str(lst.cell(lpc)),'_CC.csv']);
    end
    gdata.xaxis = xlbls;
    grandout{lpc} = gdata;
    
end % cell loop

% archive data
save (['analyzed_CC_' data '.mat'], 'grandout', 'gdata', 'out');
