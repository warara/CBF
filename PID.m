function [ u ] = PID( y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %y = squeeze(y);
    num_Agents = size(y, 1);
    num_Agents = size(y, 1);
    
    y0 = 1/(num_Agents)*ones(1, num_Agents)*y(1:num_Agents,1);

    
%     ep = (y0*ones(num_Agents,1)-y(:,1));
%     ev = -y(:,2) ;
%     P = 60;
%     D = 100;
%     u  = (P*ep+D*ev);
    u = 4*(y0-y(:,1));


end

