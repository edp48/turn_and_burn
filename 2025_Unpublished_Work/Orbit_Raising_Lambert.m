clear
clc
close all

%$$$$$$$$$$$$$$$$$$$ Set up Physical Parameters and ICs $$$$$$$$$$$$$$$$$$$

% Phsyical Constants
g=9.81;
mu = 3.9869044e14;                 % Graviational Parameter for Earth

% Spacecraft Parameters
mass_init=20000;
mf=9000; 
isp=1000;             
thrust=300;
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

p.mu = mu;
p.thrust = thrust;
p.m_dot = m_dot;
p.mass_init = mass_init;
p.tstop = tstop;
p.Isp = isp;
p.a0 = a0;
p.af = af;
p.e0 = e0;
p.ef = ef;

%$$$$$$$$$$$$$$$$$$$$ Calling the Lambert Solver $$$$$$$$$$$$$$$$$$$$$$$$$$

TOF = 19000;

[v0_lamb,vf_lamb] = Universal_Lambert(r0,rf,TOF,1,mu);

Z0 = [r0;v0_lamb];
t_int = [0 TOF];


tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
[tout,Zout] = ode113(@(t,X) Dynamic_Model(X,p),t_int,Z0,options);


dv1 = abs(norm(v0 - v0_lamb));
dv2 = abs(norm(vf - vf_lamb));

dv_tot = dv1 + dv2

dv_max = g*isp*log(mass_init/(mass_init - mf))


theta = linspace(0, 2*pi, 100);  % Angle values

% Circle 1
x1 = a0 * cos(theta);
y1 = a0 * sin(theta);

% Circle 2
x2 = af * cos(theta);
y2 = af * sin(theta);

% Plot
figure;
plot(x1, y1, 'b-', 'LineWidth', 2); hold on;
plot(x2, y2, 'r-', 'LineWidth', 2);
plot(Zout(:,1),Zout(:,2),LineWidth=3)
plot(r0(1),r0(2),'Marker','.','MarkerSize',20)
plot(rf(1),rf(2),'Marker','.','MarkerSize',20)
axis equal;
grid on;


%% Functions


function [dXdt] = Dynamic_Model(X,p)

   
    
    r = X(1:3,:);
    v = X(4:6,:);


    dXdt = zeros(size(X));

    dXdt(1,:) = X(4,:);
    dXdt(2,:) = X(5,:);
    dXdt(3,:) = X(6,:);
    

    dXdt(4:6,:) = -(p.mu./(vecnorm(r).^3)).*r;
end