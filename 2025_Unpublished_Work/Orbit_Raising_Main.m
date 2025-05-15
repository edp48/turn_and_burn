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
thrust=200;
m_dot=thrust/(g*isp);

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

% Simulation Parameters
tstop=mf/m_dot;


% Package Parameters for later use

TOF = tstop;

p.mu = mu;
p.thrust = thrust;
p.m_dot = m_dot;
p.mass_init = mass_init;
p.tstop = tstop;
p.Isp = isp;
p.a0 = a0;
p.af = af;
p.e0 = 0.1;
p.ef = 0.2;
p.t_flip = 0.437*TOF;


t_int = [0 TOF];
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,X) Open_Loop_Events(X,p));
[tout,Zout] = ode113(@(t,X) Open_Loop_Dynamics(t,X,p),t_int,X0,options);


theta = linspace(0, 2*pi, 100);  % Angle values

% Circle 1
x1 = a0 * cos(theta);
y1 = a0 * sin(theta);

% Circle 2
x2 = af * cos(theta);
y2 = af * sin(theta);

temp = 505;

% Plot
figure;
plot(x1, y1, 'b-', 'LineWidth', 2); hold on;
plot(x2, y2, 'r-', 'LineWidth', 2);
plot(Zout(:,1),Zout(:,2),LineWidth=3)
plot(r0(1),r0(2),'Marker','.','MarkerSize',20)
plot(rf(1),rf(2),'Marker','.','MarkerSize',20)
axis equal;
grid on;



if TOF>tstop
    error("TOF too long")
end


Zout_new = Zout;
Zout_new(1,:) = X0;
Zout_new(end,:) = Xf;


%$$$$$$$$$$$$$$$$$$$$$$$ Set Up and Call OptimTraj $$$$$$$$$$$$$$$$$$$$$$$$

problem.func.dynamics = @(t,X,u)(Dynamic_Model(t,X,u,p));
problem.func.pathObj = @(t,X,u)(0.0001.*t);
problem.func.bndCst = @(t0,X0,tf,Xf)(boundry_constraints(X0,Xf,p));


problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = tstop;
problem.bounds.finalTime.low = 0;
problem.bounds.finalTime.upp = tstop;

problem.bounds.control.low = -2*pi;
problem.bounds.control.upp = 2*pi;


[~, l] = min(abs(tout - p.t_flip));
one_arr = [zeros([1 l]) pi.*ones([1,length(tout)-l])];

%problem.guess.time = tout.';
%problem.guess.state = Zout_new.';
%problem.guess.control = one_arr;


problem.guess.time = tout.';
problem.guess.state = Zout.';
problem.guess.control = one_arr;

problem.options.method = 'chebyshev';
problem.options.chebyshev.nColPts = 13;
problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',500e4,...
    'tolFun',1e-10, ...
    'MaxIter',30e4);


soln = optimTraj(problem);

%save("first_output.mat","soln")


tGrid = soln(end).grid.time;
t_int = linspace(tGrid(1),tGrid(end),9000);
z = soln(end).interp.state(t_int);
u = soln(end).interp.control(t_int);
ang_int = soln(end).interp.control;

x_data = z(1,:);

y_data = z(2,:);


plot(x_data,y_data)




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

function [value,isterminal,direction] = Open_Loop_Events(X,p)

    r = X(1:3);
    r = reshape(r,1,3);




    value = norm(r) - p.af;
    isterminal = 1;
    direction = 0;




end