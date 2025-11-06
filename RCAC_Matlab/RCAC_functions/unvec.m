function [ A ] = unvec( vecA )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

l = length(vecA);
sl = sqrt(l);
A = 0;
if rem(sl,1)==0
    A = zeros(sl,sl);
    
    for ii = 1:sl
        A(:,ii) = vecA((ii-1)*sl+1 : ii*sl);
    end
    
    
    
end


end

