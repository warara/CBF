function [ dx ] = dgl_uni_four(~, x, u )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dx = [dgl_rel_uni(x(1:2), u(1)); dgl_rel_uni(x(3:4), u(2)); dgl_rel_uni(x(5:6), u(3)); dgl_rel_uni(x(7:8), u(4)) ];

end

