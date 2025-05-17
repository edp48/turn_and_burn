clear
clc
close all

%$$$$$$$$$$$$$$$$$$$ Set up Physical Parameters and ICs $$$$$$$$$$$$$$$$$$$

% Phsyical Constants
g=9.81;
mu = 3.9869044e14;                 % Graviational Parameter for Earth

% Spacecraft Parameters
mass_init=2000;
mf=1500; 
isp=1000;             


% Initial Conditions:
a0 = 8000e3;
e0 = 0.001;
I0 = 0;
omega0 = 0;
Omega0 = 0;
nu0 = 0;

[r0,v0] = orbitalelements2posvel(a0,e0,Omega0,I0,omega0,mu,nu0);
X0 = [r0;v0];

%Final Conditions
af = 42000e3;
ef = 0.001;
If = 0;
omegaf = 0;
Omegaf = 0;
nuf = 3;

[rf,vf] = orbitalelements2posvel(af,ef,Omegaf,If,omegaf,mu,nuf);
Xf = [rf;vf];


% Package Parameters for later use


p.mu = mu;
p.mass_init = mass_init;
p.Isp = isp;
p.a0 = a0;
p.af = af;
p.e0 = 0.1;
p.ef = 0.2;

% Load In Initial Guess
thrust=700;
m_dot=thrust/(g*isp);
tstop=mf/m_dot;
TOF = tstop;
p.thrust = thrust;
p.m_dot = m_dot;
p.tstop = tstop;
p.t_flip = 0.49*TOF;


t_int = [0 TOF];
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,X) Open_Loop_Events(X,p));
[tout_kmin1,Zout_kmin1] = ode113(@(t,X) Open_Loop_Dynamics(t,X,p),t_int,X0,options);

[~, l] = min(abs(tout_kmin1 - p.t_flip));
control_kmin1 = [zeros([1 l]) pi.*ones([1,length(tout_kmin1)-l])];




%$$$$$$$$$$$$$$$$$$$$$$$ Set Up and Call OptimTraj $$$$$$$$$$$$$$$$$$$$$$$$




problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = tstop;
problem.bounds.finalTime.low = 0;
problem.bounds.finalTime.upp = tstop;

problem.bounds.control.low = -2*pi;
problem.bounds.control.upp = 2*pi;

problem.options.method = 'chebyshev';
problem.options.chebyshev.nColPts = 13;
problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',500e4,...
    'tolFun',1e-10, ...
    'MaxIter',30e4);

% 67 iters at 2 newtons
% 134 iters at 1 newton

nsims = 1272;
mesh_num = 5000;

Z_Master = zeros(6,mesh_num,nsims);
Control_Master = zeros(nsims,mesh_num);
Time_Master = zeros(nsims,mesh_num);
Thrust_Master = zeros(nsims,1);
TOF_Master = zeros(nsims,1);
Lambert_Master = zeros(nsims,1);
Mdot_Master = zeros(nsims,1);

for k = 1:nsims
    
    m_dot=thrust/(g*isp);
    tstop=mf/m_dot;
    TOF = tstop;
    p.thrust = thrust;
    p.m_dot = m_dot;
    p.tstop = tstop;

    problem.func.dynamics = @(t,X,u)(Dynamic_Model(t,X,u,p));
    problem.func.pathObj = @(t,X,u)(0.0001.*t);
    problem.func.bndCst = @(t0,X0,tf,Xf)(boundry_constraints(X0,Xf,p));
    
   
    
    problem.guess.time = tout_kmin1.';
    problem.guess.state = Zout_kmin1.';
    problem.guess.control = control_kmin1;
    
    soln = optimTraj(problem);

    if soln.info.exitFlag <= 0 
        error("no convergance: " + k)
    end

    Zout_kmin1 = soln.grid.state.';
    tout_kmin1 = soln.grid.time.';
    control_kmin1 = soln.grid.control;

    Thrust_Master(k) = thrust;

    tGrid = soln(end).grid.time;
    t_int = linspace(tGrid(1),tGrid(end),mesh_num);
    z = soln(end).interp.state(t_int);
    u = soln(end).interp.control(t_int);

    Z_Master(:,:,k) = z;
    Control_Master(k,:) = u;
    Time_Master(k,:) = t_int;
    TOF_Master(k) = tout_kmin1(end);
    Mdot_Master(k) = p.m_dot;


    thrust = thrust+1;

    Lambert_Master(k)= find_lambert_soln(soln,p);



end

master.Control = Control_Master;
master.Lambert = Lambert_Master;
master.Thrust = Thrust_Master;
master.Time = Time_Master;
master.TOF = TOF_Master;
master.Z = Z_Master;
master.Mdot = Mdot_Master;

save("out3.mat","Mdot_Master","Z_Master","TOF_Master","Time_Master","Thrust_Master","Lambert_Master","Control_Master","p","master","nsims")



%%
[dv,turn] = analyze_outputs(master,p,nsims);

[Lambert_short,Lambert_long] = find_lambert_soln_alt(master,p);





%$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Plotting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
theta = linspace(0, 2*pi, 100);  % Angle values

% Circle 1
x1 = p.a0 * cos(theta);
y1 = p.a0 * sin(theta);

% Circle 2
x2 = p.af * cos(theta);
y2 = p.af * sin(theta);



% Plot
figure(1);
title("Trajectories","Interpreter","latex","FontSize",18)
plot(x1, y1, 'b-', 'LineWidth', 2); hold on;
plot(x2, y2, 'r-', 'LineWidth', 2);

axis equal;
grid on;

figure(2)
hold on
title("Control History","Interpreter","latex","FontSize",18)

for k = 1:nsims
    figure(1)
    plot(Z_Master(1,:,k),Z_Master(2,:,k))

    figure(2)
    plot(Control_Master(k,:))

end

figure(3)
plot(Thrust_Master,TOF_Master,"LineWidth",3);
title("TOF","Interpreter","latex","FontSize",18)

figure(4)

hold on
plot(Thrust_Master,dv,"DisplayName","Cont. Thrust","LineWidth",3)
plot(Thrust_Master,Lambert_short,"DisplayName","Lambert Short Way","LineWidth",3)
plot(Thrust_Master,Lambert_long,"DisplayName","Lambert Long Way","LineWidth",3)
legend
title("$\Delta V$ Values vs Thrust","Interpreter","latex","FontSize",18)

figure(5)
hold on
plot(TOF_Master,dv,"DisplayName","Cont. Thrust","LineWidth",3)
plot(TOF_Master,Lambert_short,"DisplayName","Lambert Short Way","LineWidth",3)
plot(TOF_Master,Lambert_long,"DisplayName","Lambert Long Way","LineWidth",3)
legend
title("$\Delta V$ Values vs TOF","Interpreter","latex","FontSize",18)


figure(6)
plot(Thrust_Master,turn./TOF_Master.',"LineWidth",3)
title("flip time (percent of TOF)","Interpreter","latex","FontSize",18)





%% Functions

function [dXdt] = Open_Loop_Dynamics(t,X,p)

   
    
    r = X(1:3,:);
    v = X(4:6,:);
    
    
    v_hat = v./vecnorm(v);
    

    dXdt = zeros(size(X));

    dXdt(1,:) = X(4,:);
    dXdt(2,:) = X(5,:);
    dXdt(3,:) = X(6,:);
    
    m_time = p.mass_init - p.m_dot.*t;
    
    if t<=p.t_flip
        thrust = p.thrust;
    else
        thrust = -p.thrust;
    end

    dXdt(4:6,:) = (1./m_time).*((-m_time).*((p.mu./(vecnorm(r).^3)).*r) + thrust.*v_hat);


end

function [dXdt] = Reintegrate_Dynamics(t,X,p,ang_int)

    ang = ang_int(t);
    
    r = X(1:3,:);
    v = X(4:6,:);
    
    
    v_hat = v./vecnorm(v);
    t_hat = rot_z(v_hat,ang);
    

    dXdt = zeros(size(X));

    dXdt(1,:) = X(4,:);
    dXdt(2,:) = X(5,:);
    dXdt(3,:) = X(6,:);
    
    m_time = p.mass_init - p.m_dot.*t;

    thrust = p.thrust;


    dXdt(4:6,:) = (1./m_time).*((-m_time).*((p.mu./(vecnorm(r).^3)).*r) + thrust.*t_hat);


end

function [value,isterminal,direction] = Open_Loop_Events(X,p)

    r = X(1:3);
    r = reshape(r,1,3);




    value = norm(r) - p.af;
    isterminal = 1;
    direction = 0;




end

function dv = find_lambert_soln(soln,p)

    Z = soln.grid.state.';

    r0 = Z(1,1:3);
    rf = Z(end,1:3);
    TOF = soln.grid.time(end) - soln.grid.time(1);

    [v0_lamb,vf_lamb] = Universal_Lambert(r0,rf,TOF,-1,p.mu);

    dv1 = abs(norm(Z(1,4:6) - v0_lamb));
    dv2 = abs(norm(Z(end,4:6)-vf_lamb));

    dv = dv1 + dv2;




end


function [dv,turn_time] = analyze_outputs(master,p,nsims)

% master.Control = Control_Master;
% master.Lambert = Lambert_Master;
% master.Thrust = Thrust_Master;
% master.Time = Time_Master;
% master.TOF = TOF_Master;
% master.Z = Z_Master;


    % Calculate Cont_thrust dv



    for row = 1:nsims
        dv(row) = 9.81*p.Isp.*log(p.mass_init./(p.mass_init - master.Mdot(row).*master.TOF(row)));
        idx = find(master.Control(row,:) > pi, 1);
        turn_time(row) = master.Time(row,idx);
    end






end


function [dv_short,dv_long] = find_lambert_soln_alt(master,p)


    for k = 1:length(master.Thrust)

    Z = master.Z(:,:,k)';
    Z = squeeze(Z);

    r0 = Z(1,1:3);
    rf = Z(end,1:3);
    TOF = master.TOF(k);

    [v0_short,vf_short] = Universal_Lambert(r0,rf,TOF,1,p.mu);

    dv1_short = abs(norm(Z(1,4:6) - v0_short));
    dv2_short = abs(norm(Z(end,4:6)-vf_short));
    dv_short(k) = dv1_short + dv2_short;

    [v0_long,vf_long] = Universal_Lambert(r0,rf,TOF,-1,p.mu);

    dv1_long = abs(norm(Z(1,4:6) - v0_long));
    dv2_long = abs(norm(Z(end,4:6)-vf_long));
    dv_long(k) = dv1_long + dv2_long;


    end


end



