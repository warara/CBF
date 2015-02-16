function out = CBF_calc(dgl, y0, para)

    %
    % F_r = @(v)(0.1+5*v+0.25*v^2);
    %
    % phi1= @(x,z)( 2*(x(2)-v_d)/m);
    % phi0= @(x,z)(-2*(x(2)-v_d)*F_r(x(2))/m+ep*(x(2)-v_d)^2);

    h=  @(x,z)(norm(z)-para.d);
    hf= @(x,z, dz)(h(x,z)-0.5*(dz^2)/(0.3*para.g));

    B = @(x, z, dz)(1/h(x,z)+0.5*(z'*dz/norm(z))^2);
    Bf= @(x, z, dz)(1/hf(x,z, dz));

    dBdz = @(x,z,dz)(-1/(norm(z)-d)^2*z/norm(z)+z'*(dz/norm(z)*(norm(z)-z'*z/norm(z))*(dz))/norm(z)^2);
    del_B = @(x,z, dz)([0;z'*(dz)/norm(z)^2*z; dBdz(x,z,dz)] );
    del_Bf = @(x,z,dz)(-1/hf(x,z,dz)^2*[0; -1/(0.3*para.g)*dz; z/norm(z)]);

    LfB = @(x,z, dz) (del_B (x,z, dz)'*[x(2); 0; dz]);
    LgB = @(x,z, dz) (del_B (x,z, dz)'*[0; 1/para.m; 0]);

    LfBf = @(x,z, dz) (del_Bf (x,z, dz)'*[x(2); 0; dz]);
    LgBf = @(x,z, dz) (del_Bf (x,z, dz)'*[0; 1/para.m; 0]);

    A_cbf=@(x,z, dz)[LgB(x,z, dz)];
    b_cbf=@(x,z, dz)(-LfB(x,z, dz)+gamma/B(x,z, dz));

    A_fcbf=@(x,z, dz)[LgBf(x,z, dz)];
    b_fcbf=@(x,z, dz)(-LfBf(x,z, dz)+para.gamma/Bf(x,z, dz));

    A_cc=[1; -1];
    b_cc=0.3*para.m*para.g*ones(2,1);

    H = 2*[1/para.m^2 0; 0 1e5];
    F =@(v)(-2*[F_r(v)/para.m^2; 0]);

    t_=0:para.dt:para.simTime;

    y = zeros(para.num_Agents, length(t_), 1,2);
    y(:,1,:,:) = y0;




    % y(:,1,:,1) = 100*(rand(num_Agents, 1)-0.5);
    % y(:,1,:,2) = zeros(10,1);

    options_quad = optimset('Display','off');
    options_ode45 = odeset();%('RelTol', 1e-8, 'AbsTol', 1e-8);
    h_it = zeros(length(t_), 2);
    u_it = zeros(length(t_), para.num_Agents);

    para.F_r=0;
    para.m = 1;


    H = 1;
    fminconFail = zeros(length(t_), para.num_Agents);

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

        ud = PID(y(:,i-1,:,:));
        u = zeros(para.num_Agents,1);
        for j = 1:para.num_Agents
    %         A = zeros(2*num_Agents-1,1);
    %         b = zeros(2*num_Agents-1,1);
            A = zeros(para.num_Agents+1,1);
            b = zeros(para.num_Agents+1,1);
            %h_it(i,:) = [h(y(1,i-1,:,1),(y(1,i-1,:,1)-y(3,i-1,:,1))); hf(y(1,i-1,:,1),(y(1,i-1,:,1)-y(3,i-1,:,1)),(y(1,i-1,:,2)-y(3,i-1,:,2)))];
            l = 1;         
            for k = 1:para.num_Agents
                if(j~=k)
    %                 A(l) =   A_cbf (y(j,i-1,:,:), (y(j,i-1,:,1)-y(k,i-1,:,1)), (y(j,i-1,:,2)-y(k,i-1,:,2)));
    %                 b(l) =   b_cbf (y(j,i-1,:,:), (y(j,i-1,:,1)-y(k,i-1,:,1)), (y(j,i-1,:,2)-y(k,i-1,:,2)));
                    A(l) = A_fcbf(y(j,i-1,:,:), (y(j,i-1,:,1)-y(k,i-1,:,1)), (y(j,i-1,:,2)-y(k,i-1,:,2)));
                    b(l) = b_fcbf(y(j,i-1,:,:), (y(j,i-1,:,1)-y(k,i-1,:,1)), (y(j,i-1,:,2)-y(k,i-1,:,2)));
                    l=l+1;
                end
            end
            A(end-1:end) = A_cc;
            b(end-1:end) = b_cc;
            [u(j), ~, exitflag] = fmincon(@(u)(u-ud(j))'*H*(u-ud(j)), 0, A, b, [], [], [], [], [], options_quad);

    %         zw = quadprog(H, 0, A, b-A*ud(j));
    %         u(j) = zw+ud(j);

    %         if(exitflag == -2)
    %             display('fmincon failed');
    %         end
            fminconFail(i,j) = exitflag;
            y(j,i,:,:) = simulate_step( dgl, y(j,i-1,:,:), u(j), para );
            u_it(i,:) = u;
        end
    end
    plot_states(t_(1:i-1), y(:,1:i-1,:,1), para);
    out.y = y;
    out.u = u_it;
    out.fminconFail = fminconFail;

end

function plot_states(t, y, para)

    subplot(2,1,1);
    plot(t, y)
    subplot(2,1,2);
    plot(t,  [zeros(1,length(t)); (abs(y(1,:,:,1))+abs(y(2,:,:,1)))-d;...
        (abs(y(1,:,:,1))+abs(y(2,:,:,1)))-para.d-0.5*((y(1,:,:,2)-y(2,:,:,2)).^2)/(0.3*para.g)]');

end


%     A1 = [ A_cbf(y(1,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(1,i-1,:,2));  A_cc];
%     b1 = [ b_cbf(y(1,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(1,i-1,:,2));  b_cc];
%     A2 = [ A_cbf(y(2,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(2,i-1,:,2));  A_cc];
%     b2 = [ b_cbf(y(2,i-1,:,:), norm(y(1,i-1,:,1)-y(2,i-1,:,1)), -y(2,i-1,:,2));  b_cc];

%     [u_(1),fval1,exitflag1,output1] = fmincon(@(u)(u-ud(1))'*H*(u-ud(1)), 0, A1, b1, [], [], [], [], [], options_quad);
%     [u_(2),fval2,exitflag2,output2] = fmincon(@(u)(u-ud(2))'*H*(u-ud(2)), 0, A2, b2, [], [], [], [], [], options_quad);






