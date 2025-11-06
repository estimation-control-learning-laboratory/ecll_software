function [ x_1 ] = VanDerPol( x, u, mu )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


Ts      = 0.001;

% mu      = 1;

x       = x(:);
x_1(1)  = x(1) + Ts * x(2);
x_1(2)  = x(2) + Ts * mu * (1 - x(1)^2) * x(2) - Ts * x(1) + Ts * u;

x_1     = x_1(:);

end

