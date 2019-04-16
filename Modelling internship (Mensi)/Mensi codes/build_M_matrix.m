function [M] = build_M_matrix(V,spike,t_refr,nbr_bink,bink_size,sampling_freq)

timestep = 1e3/sampling_freq;
M = zeros(nbr_bink,length(V)+200000);
spike = spike + round((t_refr)/timestep);

for i=1:length(spike)
    for j=1:nbr_bink
        if(j==1)
            M(j,spike(i):spike(i)+bink_size(j)-1) = ...
                M(j,spike(i):spike(i)+bink_size(j)-1) + ones(1,bink_size(j));
        else
            M(j,spike(i)+sum(bink_size(1:(j-1))):spike(i)+sum(bink_size(1:(j-1)))+bink_size(j)-1) = ...
                M(j,spike(i)+sum(bink_size(1:(j-1))):spike(i)+sum(bink_size(1:(j-1)))+bink_size(j)-1) + ones(1,bink_size(j));
        end
    end
end

M = M(:,1:length(V));





