function [ dx ] = dgl_rel_uni(~, x, u, param )
%DGL_REL_UNI Summary of this function goes here
%   Detailed explanation goes here

    dx = zeros(5,1);
    dx(1) = x(4)*cos(x(3));
    dx(2) = x(4)*sin(x(3));
    dx(3) = x(5);
    dx(4) = -param.F_r/param.m+u(1)/param.m;
    dx(5) = u(2);

end

