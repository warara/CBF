para.g = 9.81;
para.m = 1;
para.gamma = 10;
para.ep = 10;
para.d = 1;
para.dt = 1e-2;
para.simTime = 10;

para.num_Agents = 4;
y0 = zeros(para.num_Agents, 2);
y0(1,:) = [ 9; 5];
y0(3,:) = [ 3; -1];
y0(2,:) = [-9; 1e-3];
y0(4,:) = [22; 0];
% y0(5,:) = [-33; -3];
% y0(6,:) = [-22; -7];
% y0(7,:) = [-15; -3];
% y0(8,:) = [ 25; 10];

out = CBF_calc(@dgl_uni, y0, para);