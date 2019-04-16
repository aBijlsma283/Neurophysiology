dt1 = 0.05/1000;
T = 2;
Nstep1 = round(T/dt1);

spiketrain1 = zeros(1,Nstep1);
Nspike = 50;

spikeindices = randi(length(spiketrain1), [1,Nspike]);
spiketrain1(spikeindices) = 1;

dt2 = 0.5/1000;
factor = round(dt2/dt1);
Nstep2= round(T/dt2);

spiketrain2 = nan(1,Nstep2);

for n = 1:Nstep2
    spiketrain2(n) = sum(spiketrain1((n-1)*factor+1:factor*n));
end

figure
plot(1:Nstep1, spiketrain1, '.')
hold all
plot(factor*(1:Nstep2)-factor/2, 2*spiketrain2, '.')
ylim([0.5 2.5])