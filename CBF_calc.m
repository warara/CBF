g = 9.81;
m = 1;
v0 = 13.89;
v_d = 24;
v_d2 = v0;
gamma = 1;
ep = 10;
d = 1;
para.dt = 1e-2;
num_Agents = 3;


%
% F_r = @(v)(0.1+5*v+0.25*v^2);
%
% phi1= @(x,z)( 2*(x(2)-v_d)/m);
% phi0= @(x,z)(-2*(x(2)-v_d)*F_r(x(2))/m+ep*(x(2)-v_d)^2);

h=  @(x,z)(z-d);
hf= @(x,z, dz)(z-d-0.5*(dz)^2/(0.3*para.m));
B = @(x,z, dz)(1/h(x,z)+0.5*(dz)^2);
del_B =  @(x,z, dz)([0;x(2); -1/(z-d)^2] );

%Bf = @(x,z,dz)();

LfB = @(x,z, dz) (del_B (x,z, dz)'*[x(2); -0/m; dz]);
LgB = @(x,z, dz) (del_B (x,z, dz)'*[0; 1/m; 0]);

A_cbf=@(x,z, dz)[LgB(x,z, dz)];
b_cbf=@(x,z, dz)(-LfB(x,z, dz)+gamma/B(x,z, dz));

A_cc=[1; -1];
b_cc=0.3*m*g*ones(2,1);


H = 2*[1/m^2 0; 0 1e5];
F =@(v)(-2*[F_r(v)/m^2; 0]);

t_=0:para.dt:10;

y = zeros(num_Agents, length(t_), 1,2);
y(1,1,:,:) = [ 9; 0];
y(2,1,:,:) = [ 1; -3];
y(3,1,:,:) = [-9; 0];

dgl = @dgl_rel_uni;
options_quad = optimset('Display','off');
options_ode45 = odeset();%('RelTol', 1e-8, 'AbsTol', 1e-8);
h_it = zeros(length(t_), 2);
u_it = zeros(length(t_), 2);

para.F_r=0;
para.m = 1;


H = 1;
fminconFail = 0;

% H = 2*[1/para.m^2 0; 0 1e5];
% F =@(v)(-2*[0; 0]);
%
% dVdx1 = @(x, z)(norm(z)-d)*x(2)/norm(z)+2*x(1)*x(2)/norm(z)*((x(2)*norm(z)-x(1)/norm(z))/norm(z)^2);
% dVdx2 = @(x, z)(2*x(1)^2*x(2)/norm(z)^2);

% dV = @(x,z)[dVd1(x,z); dVdx2(x,z)];
%
% f =@(x) [x(2); -F_r(x(2))/para.m];
% g = [0; para.m];

for i=2:length(t_)

    
    dy =  0;
    ud = PID(y(:,i-1,:,:));
   
    u = zeros(num_Agents,1);
    for j = 1:num_Agents
        A = zeros(num_Agents-1,1);
        b = zeros(num_Agents-1,1);
        l = 1;
        for k = 1:num_Agents
            if(j~=k)
                A(l) = A_cbf(y(j,i-1,:,:), norm(y(j,i-1,:,1)-y(k,i-1,:,1)), -y(j,i-1,:,2));
                b(l) = b_cbf(y(j,i-1,:,:), norm(y(j,i-1,:,1)-y(k,i-1,:,1)), -y(j,i-1,:,2));
                l=l+1;
            end
        end
        A(end-1:end) = A_cc;
        b(end-1:end) = b_cc;
        [u(j), ~, exitflag] = fmincon(@(u)(u-ud(j))'*H*(u-ud(j)), 0, A, b, [], [], [], [], [], options_quad);
    end

    y(:,i,:,:) = simulate_step( dgl, y(:,i-1,:,:), u, para );
        
end
hold on
subplot(2,1,1);

plot(t_(1:i-1), [y(:,1:i-1,:,1);])
subplot(2,1,2);
plot(t_(1:i-1),  [zeros(1,i-1); abs(y(1,1:i-1,:,1)- y(2,1:i-1,:,1))-d; abs(y(3,1:i-1,:,1)- y(2,1:i-1,:,1))-d]');





%     A1 = [ A_cbf(y(1,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(1,i-1,:,2));  A_cc];
%     b1 = [ b_cbf(y(1,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(1,i-1,:,2));  b_cc];
%     A2 = [ A_cbf(y(2,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(2,i-1,:,2));  A_cc];
%     b2 = [ b_cbf(y(2,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(2,i-1,:,2));  b_cc];

%      [u_(1),fval1,exitflag1,output1] = fmincon(@(u)(u-ud(1))'*H*(u-ud(1)), 0, A1, b1, [], [], [], [], [], options_quad);
%      [u_(2),fval2,exitflag2,output2] = fmincon(@(u)(u-ud(2))'*H*(u-ud(2)), 0, A2, b2, [], [], [], [], [], options_quad);






