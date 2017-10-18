function [ a ] = compute_homography( a )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    H = [1     0.1   0;
     0.1   1     0;
     0.004 0.002 1 ];
    a = H*a;
end

