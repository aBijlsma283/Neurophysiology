function [theta nbr_iter hist_L hist_theta] = compute_theta(X,X_spike,sum_X_spike,theta_0,tol_theta,max_iter)

X = X'; X_spike = X_spike'; theta_0 = theta_0';

X_spike_theta = X_spike*theta_0;
X_theta = X*theta_0; 
expXtheta = exp(X_theta);

L = sum(X_spike_theta) - sum(expXtheta);
G = sum_X_spike' - (X'*expXtheta);
temp = nan(size(X'));
for j=1:size(X,2)
    temp(j,:) = X(:,j).*expXtheta;
end
H = -(temp*X);

hist_theta = nan(max_iter,length(theta_0));
hist_L = nan(max_iter,1);
hist_theta(1,:) = theta_0;
hist_L(1) = L;
theta = theta_0;

for i=1:max_iter
    if(i>3 && abs((L-hist_L(i-2))/hist_L(i-2)) < tol_theta)
        nbr_iter = i;
        hist_theta(i,:) = theta;
        hist_L(i) = L;
        return;
    else
        nbr_iter = i;
        theta = theta - ((inv(H))*G);
        X_spike_theta = X_spike*theta;
        X_theta = X*theta;
        expXtheta = exp(X_theta);
        
        L = sum(X_spike_theta) - sum(expXtheta);
        G = sum_X_spike' - (X'*expXtheta);
        temp = nan(size(X'));
        for j=1:size(X,2)
            temp(j,:) = X(:,j).*expXtheta;
        end
        H = -(temp*X);
        
        hist_theta(i,:) = theta;
        hist_L(i) = L;
    end
   
end
