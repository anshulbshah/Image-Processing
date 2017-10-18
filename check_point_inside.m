function [ b ] = check_point_inside( a )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if(a(1) >= 1 && a(1) <= 800 && a(2) >= 1 && a(2) <= 600)
        b = 1
    else
        b = 0
    end
end

