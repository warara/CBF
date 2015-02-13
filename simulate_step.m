function [ y_new ] = simulate_step( dgl, y, u, para )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    y=squeeze(y);
    num_Agents = size(y, 1);
    y_new = zeros(num_Agents,2);
    for i=1:num_Agents
        [~, Y] = ode45( @(t,x)(dgl(t, x, u(i), para)), [0 para.dt], y(i,:));
        y_new(i,:) = Y(end,:);
    end


end

