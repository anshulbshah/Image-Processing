function [ mu,N ] = mu_2( t,im )
%MU_1 Summary of this function goes here
%   Detailed explanation goes here
sum = 0;
for j = t+1:255
    sum = sum + j*size(find(im==j),1);
end
mu = sum/size(find(im>t),1);
N = size(find(im>t),1);
