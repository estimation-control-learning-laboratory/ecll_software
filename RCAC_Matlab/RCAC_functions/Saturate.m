function [ usat ] = Saturate( u, llim, ulim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

usat = max(u, llim);
usat = min(usat, ulim);


end

