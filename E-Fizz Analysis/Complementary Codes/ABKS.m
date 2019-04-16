function [h,p,ksstat,cv] = ABKS(data)

data(isnan(data)) = [];
Ndata = (data-mean(data))/std(data);
[h,p,ksstat,cv] = kstest(Ndata);