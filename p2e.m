function [ proj ] = p2e( proj )
%   Projective to Eucledian
    if(size(proj,1) == 3)
        proj = [proj(1,:)./proj(3,:) ; proj(2,:)./proj(3,:)];
    else
        proj = [proj(1,:)./proj(4,:) ; proj(2,:)./proj(4,:) ; proj(3,:)./proj(4,:)];
end

