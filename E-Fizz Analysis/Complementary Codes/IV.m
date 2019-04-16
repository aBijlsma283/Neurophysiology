%% IV 
%% Made by Ate Bijlsma, s4212215, ate.bijlsma@student.ru.nl
% makes an IV curve from the grouped data. As the holding potentials of
% each experiment differ we need to align the different potentials with
% each other. For now this is done by hand.

%% Aligning
close all

gdata.IVwith(isnan(gdata.IVwith)) = 0;
test = gdata.IVwith;

plus = nan(15,1);
newtest = [test(:,1) plus test(:,2) plus test(:,3) plus test(:,4) plus test(:,5) plus test(:,6) plus test(:,7) plus test(:,8) plus test(:,9) plus test(:,10) plus test(:,11) plus test(:,12)];
newtest = [plus plus plus plus plus newtest plus plus plus plus plus plus plus plus plus plus];

newtest(2,:) = [nan nan nan nan nan newtest(2,1:28) nan nan nan nan nan];
newtest(6,:) = [nan nan nan nan nan newtest(6,1:28) nan nan nan nan nan];
newtest(7,:) = [nan nan nan nan nan newtest(7,1:28) nan nan nan nan nan];
newtest(8,:) = [nan nan nan nan nan newtest(8,1:28) nan nan nan nan nan];
newtest(9,:) = [nan nan nan nan nan newtest(9,1:28) nan nan nan nan nan];
newtest(14,:) = [nan nan nan nan nan newtest(14,1:28) nan nan nan nan nan];

newtest(4,:) = [newtest(4,6:28) nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan];
newtest(5,:) = [newtest(5,6:28) nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan];

newtest(10,:) = [nan nan nan nan nan nan nan nan newtest(10,6:28) nan nan nan nan nan nan nan];

newtest(11,:) = [nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan newtest(11,6:28)];
newtest(15,:) = [nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan newtest(15,6:28)];

IVwith = newtest;

IVstim = -63:1:-26;

gdata.IVwithout(isnan(gdata.IVwithout)) = 0;
test = gdata.IVwithout;

plus = nan(15,1);
newtest = [test(:,1) plus test(:,2) plus test(:,3) plus test(:,4) plus test(:,5) plus test(:,6) plus test(:,7) plus test(:,8) plus test(:,9) plus test(:,10) plus test(:,11) plus test(:,12)];
newtest = [plus plus plus plus plus newtest plus plus plus plus plus plus plus plus plus plus];

newtest(2,:) = [nan nan nan nan nan newtest(2,1:28) nan nan nan nan nan];
newtest(6,:) = [nan nan nan nan nan newtest(6,1:28) nan nan nan nan nan];
newtest(7,:) = [nan nan nan nan nan newtest(7,1:28) nan nan nan nan nan];
newtest(8,:) = [nan nan nan nan nan newtest(8,1:28) nan nan nan nan nan];
newtest(9,:) = [nan nan nan nan nan newtest(9,1:28) nan nan nan nan nan];
newtest(14,:) = [nan nan nan nan nan newtest(14,1:28) nan nan nan nan nan];

newtest(4,:) = [newtest(4,6:28) nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan];
newtest(5,:) = [newtest(5,6:28) nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan];

newtest(10,:) = [nan nan nan nan nan nan nan nan newtest(10,6:28) nan nan nan nan nan nan nan];

newtest(11,:) = [nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan newtest(11,6:28)];
newtest(15,:) = [nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan newtest(15,6:28)];

IVwithout = newtest;

%% Cakculation of the nan means.
for ctr = 1:38
    stdwith(1,ctr) = nanstd(IVwith(:,ctr))/(nnz(~isnan(IVwith(:,ctr))));
    stdwithout(1,ctr) = nanstd(IVwithout(:,ctr))/(nnz(~isnan(IVwithout(:,ctr))));
end

meanIVwith = nanmean(IVwith,1);
meanIVwithout = nanmean(IVwithout,1);

IVstim(isnan(meanIVwith)) = [];
meanIVwithout(isnan(meanIVwith)) = [];
meanIVwith(isnan(meanIVwith)) = [];
% 

%% Main Plot

subplot(4,4,1)

difconditions = abs(meanIVwithout-meanIVwith);

hold on
plot(IVstim,-meanIVwith,'r.')
plot(IVstim,-meanIVwithout,'k.');
plot(IVstim,difconditions,'b.');
  
f = fit(IVstim',-meanIVwith','smoothingspline','SmoothingParam',0.1);
r1 = plot(f,'r');
hold on

f = fit(IVstim',-meanIVwithout','smoothingspline','SmoothingParam',0.1);
r2 = plot(f,'k');
hold on

f = fit(IVstim',difconditions','smoothingspline','SmoothingParam',0.1);
r3 = plot(f,'b');
legend off

[h,p] = ttest(-meanIVwith,-meanIVwithout);

ylabel('Current [nA]','FontSize',8)
title(['Average I/V curve VC-step. p = ' num2str(p)],'FontSize',12)
xlim([IVstim(1) IVstim(end)])
legend([r1 r2 r3],{'Magnet ON','Magnet OFF','Residual'})

IVstim = -63:1:-26;

for ctr = 1:size(IVwith,1)
    subplot(4,4,ctr+1)
    hold on
    difference(ctr,:) = abs(gdata.IVwithout(ctr,:)-gdata.IVwith(ctr,:));
    
    plot(gdata.stimwith(ctr,:),-gdata.IVwith(ctr,:),'r.')
    plot(gdata.stimwith(ctr,:),-gdata.IVwithout(ctr,:),'k.');
    plot(gdata.stimwith(ctr,:),difference(ctr,:),'b.');
    
    f = fit(gdata.stimwith(ctr,:)',-gdata.IVwith(ctr,:)','smoothingspline','SmoothingParam',0.1);
    r1 = plot(f,'r');
    f = fit(gdata.stimwith(ctr,:)',-gdata.IVwithout(ctr,:)','smoothingspline','SmoothingParam',0.1);
    r2 = plot(f,'k');
    f = fit(gdata.stimwith(ctr,:)',difference(ctr,:)','smoothingspline','SmoothingParam',0.1);
    r3 = plot(f,'b');
    
    [h(ctr),p(ctr)] = ttest(difference(ctr,:),zeros(1,20));
    
    legend off
    xlabel('Hold potential [mV]','FontSize',8)
    ylabel('Current [nA]','FontSize',8)
    title(['Exp ' num2str(ctr) ' I/V curve. p = ' num2str(p(ctr))],'FontSize',12)
    xlim([gdata.stimwith(ctr,1) gdata.stimwith(ctr,end)])
end

save('IV_Data', 'gdata', 'IVwith', 'IVwithout')
set(gcf,'Position', get(0, 'Screensize')); 
savefig('IV Curve')
  