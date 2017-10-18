function [ proj ] = p2e( proj )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    if(size(proj,1) == 3)
        proj = [proj(1,:)./proj(3,:) ; proj(2,:)./proj(3,:)];
    else
        proj = [proj(1,:)./proj(4,:) ; proj(2,:)./proj(4,:) ; proj(3,:)./proj(4,:)];
end

