function [CoFac1, CoFac2, CoFac] = calccofac_symmetric_ignoredoublespikes(eventtimes1, eventtimes2, rate1, rate2, precision)

[~, ~, ~, CoFac1] = loop_over_spikes(eventtimes1, eventtimes2, rate2, precision);
[~, ~, ~, CoFac2] = loop_over_spikes(eventtimes2, eventtimes1, rate1, precision);

CoFac = nanmean([CoFac1, CoFac2]);


    function [Ncoinc, Nbar, Ncurl, CoFac] = loop_over_spikes(eventtimes1, eventtimes2, rate2, precision)
        N1=length(eventtimes1);
        N2=length(eventtimes2);
        
        coinc=[];
        for n=1:N1
            clear a
            a=find(eventtimes2>=eventtimes1(n)-precision&eventtimes2<=eventtimes1(n)+precision);
            if isempty(a)
            elseif length(a)==1
                coinc=[coinc eventtimes1(n)];
            elseif length(a)>1
                % disp('two coincident spikes - ignore')
                % continue
                coinc=[coinc eventtimes1(n)];
            end
        end
        
        Ncoinc=length(coinc);
        Nbar=2*rate2*precision*N1;
        Ncurl=1-2*rate2*precision;
        CoFac=(Ncoinc-Nbar)/((1/2)*(N1+N2))*(1/Ncurl);
        %         if (CoFac-1)>10*eps
        if (CoFac-1)>.001
            disp('Gamma>1, replace with NaN')
            CoFac = NaN;
        end
    end

end