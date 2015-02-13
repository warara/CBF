function [ output_args ] = iif( condition, a, b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if(condition)
        output_args = a;
    else
        output_args = b;
    end

end

