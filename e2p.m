function [ euc ] = e2p( euc )
%   Eucledian to projective
    euc = [euc; ones(size(euc,2),1)'];

end

