function [ u ] = PID( y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %y = squeeze(y);
    num_Agents = size(y, 1);
    if(nargin<2)
        num_Agents = size(y, 1);
        y0 = 1/(num_Agents+1)*ones(1, num_Agents)*y(1:num_Agents,1,:,1);
    end
    
    ep = (y0*ones(num_Agents,1)-y(:,1));
    ev = -y(:,2) ;
    P = 80;
    D = 120;
    u  = (P*ep+D*ev);



end
