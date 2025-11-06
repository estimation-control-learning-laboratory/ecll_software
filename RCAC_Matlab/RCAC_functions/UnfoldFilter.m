function [ N, Nsum ] = UnfoldFilter( Nu, lu )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lz  = size(Nu, 1);
nf  = size(Nu, 2)/lu;

N   = zeros(lz, lu, nf);
Nsum= zeros(lz, lu);
for ii = 1:nf
    N(:,:,ii) = Nu(:, (ii-1)*lu+1 : ii*lu);
    
    Nsum = Nsum + Nu(:, (ii-1)*lu+1 : ii*lu);
end



end

