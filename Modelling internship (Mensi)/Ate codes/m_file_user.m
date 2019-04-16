for paramctr = 1:4
    [vv_mp(paramctr,:),Mc(paramctr,1),spikecount_mp(paramctr,1)] = IF_eta_modified(m_I(paramctr,:),m_param(paramctr,:),m_eta(paramctr,:)',m_gamma(paramctr,:)',m_sampling_freq(paramctr,:),m_Spikecount(paramctr,:));
end