function [ dx ] = dgl_rel_uni(~, x, u, param )
%DGL_REL_UNI Summary of this function goes here
%   Detailed explanation goes here

    dx = zeros(1,1);
    dx(1) = u;
    
end

