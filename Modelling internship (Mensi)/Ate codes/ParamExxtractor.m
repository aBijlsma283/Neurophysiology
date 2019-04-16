warning off

analyzedfiles = dir ('*.mat');
% param = zeros(size(analyzedfiles,1),7);
n = 0;
for paramctr = 1:size(analyzedfiles,1)
    load(analyzedfiles(paramctr).name(1:end-4))
    if length(I) == 7200000
        m_I(paramctr,:) = I;
    elseif length(I) == 3600000
        m_I(paramctr,:) = [I;I];
    end
    m_param(paramctr,:) = [C g_l E_l E_reset t_refr v0 DeltaV];
    m_eta(paramctr,:) = eta;
    m_gamma(paramctr,:) = gamma;
    m_sampling_freq(paramctr,:) = sampling_freq;
    m_Spikecount(paramctr,:) = nbr_spikes;
    if length(I) == 7200000
        m_hidden_state(paramctr,:) = hidden_state;
    elseif length(I) == 3600000
        m_hidden_state(paramctr,:) = [hidden_state;hidden_state];
    end
    m_settings{paramctr,1} = settings;
    r_v(paramctr,:) =  cell_response;
    n = n + 1;
    disp(n)
end
keyboard

