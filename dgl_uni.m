function [ dx ] = dgl_rel_uni(~, x, u, param )
%DGL_REL_UNI Summary of this function goes here
%   Detailed explanation goes here

    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -param.F_r/param.m+u/param.m;

end

