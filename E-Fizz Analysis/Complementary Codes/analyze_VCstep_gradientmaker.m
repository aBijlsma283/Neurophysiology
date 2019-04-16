function [gradient]=analyze_VCstep_gradientmaker(color1,color2,count)
%% Made by Ate Bijlsma, s4212215, ate.bijlsma@student.ru.nl

% Anlyze_VCstep_gradientmaker allows you to make a gradient between 2
% colors. 
% Input is the color where the gradient has to start (color1) where the
% gradient has to end (color2) and the amount of points in the gradient
% (count).

% Colors have to be given in Intensity values, not RGB format!

% Calculating stepsize
r=(color2(1)-color1(1))/(count-1);
g=(color2(2)-color1(2))/(count-1);
b=(color2(3)-color1(3))/(count-1);

% prelocation of the gradient
gradient=zeros(count,3);

%for each color step, increase/reduce the value of Intensity data.
for ctr=1:count
    gradient(ctr,1) = color1(1)+r*(ctr-1);
    gradient(ctr,2) = color1(2)+g*(ctr-1);
    gradient(ctr,3) = color1(3)+b*(ctr-1);
end