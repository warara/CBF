para.g = 9.81;
para.m = 1;
para.gamma = 1;
para.ep = 10;
para.d = 1;
para.dt = 1e-2;
para.simTime = 5;

para.num_Agents = 3;
y0 = zeros(para.num_Agents, 1);
y0(1,:) = [ 9;];
y0(2,:) = [ 1;];
y0(3,:) = [-9;];
%y0(4,:) = [22; 0];
% y0(5,:) = [-33; -3];
% y0(6,:) = [-22; -7];
% y0(7,:) = [-15; -3];
% y0(8,:) = [ 25; 10];

out = CBF_calc(@dgl_uni_vel, y0, para);