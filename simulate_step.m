function [ y_new ] = simulate_step( dgl, y, u, para )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [~, Y] = ode45( @(t,x)(dgl(t, x, u, para)), [0 para.dt], squeeze(y));
    y_new = Y(end,:);



end

