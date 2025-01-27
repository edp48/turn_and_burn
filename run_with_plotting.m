%% Developed by Erik Payton and Isaac Beaton F23 For use with "Isaac_Model_Opt.m"
clear
clc
close all
 
tic

%load results.mat
%IMPORTANT NOTE: If an error occurs in the expected trajectory with a sharp
%change and a lot of oscilation in the thrust, check the velocity at that
%point, as if the velocity vector gains a magnitude of zero, then the
%simulation will break, you have already spent two seperate times
%discovering this!!



% ICS FOR TRAILING
% p.t_flip=0.53*tstop;
% th=deg2rad(135);
% 368N

%% set up ET (sec after J2000)
et0 = cspice_str2et( '2024 FEB 13' );
ets = (0:1:28*86400)+et0;
[state,lt] = cspice_spkezr( 'MOON',ets,'J2000','NONE' , 'EARTH' );

moon_dat.moon_x=1000.*state(1,:).';
moon_dat.moon_y=1000.*state(2,:).';
moon_dat.moon_z=0*1000.*state(3,:).';
moon_dat.moon_vx=1000.*state(4,:).';
moon_dat.moon_vy=1000.*state(5,:).';
moon_dat.moon_vz=0*1000.*state(6,:).';
moon_dat.time = (0:1:28*86400);

r_em=[moon_dat.moon_x moon_dat.moon_y moon_dat.moon_z].'; %Position of the moon wrtearth in m


funs.Intx = griddedInterpolant(moon_dat.time,moon_dat.moon_x);
funs.Inty = griddedInterpolant(moon_dat.time,moon_dat.moon_y);
funs.Intz = griddedInterpolant(moon_dat.time,moon_dat.moon_z);

funs.Intvx = griddedInterpolant(moon_dat.time,moon_dat.moon_vx);
funs.Intvy = griddedInterpolant(moon_dat.time,moon_dat.moon_vy);
funs.Intvz = griddedInterpolant(moon_dat.time,moon_dat.moon_vz);







%%


% Universal Constants
g=9.81;
MuE = 3.9869044e14;                 % Graviational Parameter for Earth in useful units 
MuM = 4.9048695e12;                 % Gravitational parameter for the Moon
MuS = 1.32712440018e20;             % Graviational parameter for the Sun

%Spacecraft Parameters
a0=0.0017*g;
mass_init=20000;
mf=9000; %mass fuel
isp=1000;                      % m/s^2 initial accel of the spacecraft...used to find mdot
thrust=300;%mass_init*a0;              % compute the constant thrust given by the initial mass/accel of the spacecraft
m_dot=thrust/(g*isp);


%Orbital Parameters
Re=6378176;                         % Radius of the earth [m]
Rm=1740000;                         % Radius of the moon 


h=35000000;% GEO orbit [m]
h_lunar = 100000;                   % LLO orbit [m]
t_period=2*pi*sqrt(((Re+h)^3)/MuE);
    % formula for velocity of a circular orbit at this r0

%Final Conditions
rf=Rm+h_lunar;
vf=(sqrt(MuM/(Rm+h_lunar))); %final velocity wrt the moon

%Simulation Parameters

tstop=mf/m_dot;

t_flip=0.515*tstop;
th=deg2rad(112);

guess_init = [t_flip th];


%% Package the needed Paramaters

p.Rm = Rm;
p.h_lunar = h_lunar;
p.Re = Re;
p.h = h;
p.MuE = MuE;
p.MuM = MuM;
p.thrust = thrust;
p.m_dot = m_dot;
p.mass_init = mass_init;
p.tstop = tstop;
p.Isp = isp;


%% First Run of the Optimizer to get t_vec output

thrust_v(1) = thrust;
options = optimset('MaxFunEval',1000,'MaxIter',20,'Display','iter','TolFun',1e-15,'TolX',1e-15);

[opt_initial_vec(:,1),fval,exitflag,output]= fminsearch(@(guess_init) cost_function_matlab(guess_init,p,funs),guess_init,options);

[TOF_v(1),m_fuel_v(1),dV_v(1),flip_t_v(1),e_s_v(1),leadvtrail(1)] = extract_info(opt_initial_vec(:,1),p,funs);


%% Run with small purterbation
n=100;

for j=2:50

% Universal Constants
g=9.81;
MuE = 3.9869044e14;                 % Graviational Parameter for Earth in useful units 
MuM = 4.9048695e12;                 % Gravitational parameter for the Moon
MuS = 1.32712440018e20;              % Graviational parameter for the Sun

%Spacecraft Parameters

mass_init=20000;
mf=9000; %mass fuel
isp=1000;                      % m/s^2 initial accel of the spacecraft...used to find mdot
thrust=thrust+2;              % compute the constant thrust given by the initial mass/accel of the spacecraft
thrust_v(j) = thrust;
m_dot=thrust/(g*isp);
p.thrust = thrust;
p.m_dot = m_dot;


%Orbital Parameters
Re=6378176;                         % Radius of the earth [m]
Rm=1740000;                         % Radius of the moon 


h=35000000;% GEO orbit [m]
h_lunar = 100000;                   % LLO orbit [m]
t_period=2*pi*sqrt(((Re+h)^3)/MuE);
    % formula for velocity of a circular orbit at this r0

%Final Conditions
rf=Rm+h_lunar;
vf=(sqrt(MuM/(Rm+h_lunar))); %final velocity wrt the moon

%Simulation Parameters

p.tstop=mf/m_dot;

j
options = optimset('MaxFunEval',5000,'MaxIter',200,'Display','iter','TolFun',1e-8,'TolX',1e-8);

[opt_initial_vec(:,j),fval,exitflag,output]= fminsearch(@(guess_init) cost_function_matlab(guess_init,p,funs),guess_init,options);

guess_init = opt_initial_vec(:,j);

[TOF_v(j),m_fuel_v(j),dV_v(j),flip_t_v(j),e_s_v(j),leadvtrail(j)] = extract_info(opt_initial_vec(:,j),p,funs);

% if leadvtrail ==1 
%     error("ghueg")
% end



end

%%

t_flip=0.533*tstop;
th=2.318;
thrust = 349;

guess_init = [t_flip th]

for j=51:100

% Universal Constants
g=9.81;
MuE = 3.9869044e14;                 % Graviational Parameter for Earth in useful units 
MuM = 4.9048695e12;                 % Gravitational parameter for the Moon
MuS = 1.32712440018e20;              % Graviational parameter for the Sun

%Spacecraft Parameters

mass_init=20000;
mf=9000; %mass fuel
isp=1000;                      % m/s^2 initial accel of the spacecraft...used to find mdot
thrust=thrust+1;              % compute the constant thrust given by the initial mass/accel of the spacecraft
thrust_v(j) = thrust;
m_dot=thrust/(g*isp)
p.thrust = thrust;
p.m_dot = m_dot;


%Orbital Parameters
Re=6378176;                         % Radius of the earth [m]
Rm=1740000;                         % Radius of the moon 


h=35000000;% GEO orbit [m]
h_lunar = 100000;                   % LLO orbit [m]
t_period=2*pi*sqrt(((Re+h)^3)/MuE);
    % formula for velocity of a circular orbit at this r0

%Final Conditions
rf=Rm+h_lunar;
vf=(sqrt(MuM/(Rm+h_lunar))); %final velocity wrt the moon

%Simulation Parameters

p.tstop=mf/m_dot;


j
options = optimset('MaxFunEval',5000,'MaxIter',200,'Display','iter','TolFun',1e-8,'TolX',1e-8);

[opt_initial_vec(:,j),fval,exitflag,output]= fminsearch(@(guess_init) cost_function_matlab(guess_init,p,funs),guess_init,options);

guess_init = opt_initial_vec(:,j);

[TOF_v(j),m_fuel_v(j),dV_v(j),flip_t_v(j),e_s_v(j),leadvtrail(j)] = extract_info(opt_initial_vec(:,j),p,funs);

if leadvtrail ==1 
    error("ghueg")
end



end

save("runboi2.mat")

%% Plotting the Data
TOF_days = TOF_v/86400;

figure(1)
scatter(thrust_v,TOF_days,"filled")
title("Time of Flight vs Thrust",'Interpreter','latex','FontSize',20)
xlabel("Thrust [N]",'Interpreter','latex','FontSize',18)
ylabel("Time of Flight [days]", 'Interpreter','latex',FontSize=18)

figure(2)
scatter(thrust_v,m_fuel_v,"filled")
title("Prop Mass Burned vs Thrust",'Interpreter','latex','FontSize',20)
xlabel("Thrust [N]",'Interpreter','latex','FontSize',18)
ylabel("Prop Mass [kg]", 'Interpreter','latex',FontSize=18)

figure(3)
scatter(thrust_v,dV_v,"filled")
title("$\Delta V$ vs Thrust",'Interpreter','latex','FontSize',20)
xlabel("Thrust [N]",'Interpreter','latex','FontSize',18)
ylabel("$\Delta V$ [$\frac{m}{s}$]", 'Interpreter','latex',FontSize=18)

flip_perc = flip_t_v./TOF_v;

figure(4)
scatter(thrust_v,flip_perc,"filled")
title("Flip Time vs Thrust",'Interpreter','latex','FontSize',20)
xlabel("Thrust [N]",'Interpreter','latex','FontSize',18)
ylabel("Flip Time [/% of TOF]", 'Interpreter','latex',FontSize=18)

figure(5)
scatter(thrust_v,e_s_v,"filled")
title("Spacecraft Eccentricity vs Thrust",'Interpreter','latex','FontSize',20)
xlabel("Thrust [N]",'Interpreter','latex','FontSize',18)
ylabel("Eccentricity", 'Interpreter','latex',FontSize=18)





%% Functions

function [TOF,m_fuel,dV,flip_t,e_s,leadvtrail] = extract_info(opt_inputs,p,funs)

    % Run the Simulation
    g=9.81;
    MuM = 4.9048695e12;                 % Gravitational parameter for the Moon
    MuE = 3.9869044e14; 
    %Orbital Parameters
    Re=6378176;                         % Radius of the earth [m]
    Rm=1740000;                         % Radius of the moon 
    %h=400000; %ALSO CHANGE IN COST
    h=35000000;% LEO orbit [m]
    h_lunar = 100000;                   % LLO orbit [m]
    rf=Rm+h_lunar;
    p.t_flip=opt_inputs(1);
    angle_set=opt_inputs(2);
    r0=(Re+h)*[cos(angle_set) sin(angle_set) 0].';
    v0=(sqrt(MuE/(Re+h)))*[-sin(angle_set) cos(angle_set) 0].';



    X0 = [r0;v0];
    t_span = [0 p.tstop];
    
    tol = 1e-11;
    
    options = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,X) turn_and_burn_events(t,X,funs,p));
    
    [tout,Xout] = ode113(@(t,X) turn_and_burn_dynamics(t,X,p,funs),t_span,X0,options);

    
    r_s_e = Xout(:,1:3);
    v_s_e = Xout(:,4:6);


%Interpolate to get the moon's geo-centric pos/vel at the correct time
%points


r_m_e(:,1) = funs.Intx(tout);
r_m_e(:,2) = funs.Inty(tout);
r_m_e(:,3) = funs.Intz(tout);

v_m_e(:,1) = funs.Intvx(tout);
v_m_e(:,2) = funs.Intvy(tout);
v_m_e(:,3) = funs.Intvz(tout);

% Calculate the spacecraft's lunar-centric pos/vel
r_s_m = r_s_e-r_m_e;
v_s_m = v_s_e-v_m_e;


%Calculate eccentricity about the moon
h=cross(r_s_m(end,:),v_s_m(end,:));

e_vec=(cross(v_s_m(end,:),h))/p.MuM - r_s_m(end,:)/norm(r_s_m(end,:));

e_s=norm(e_vec); %e of the spacecraft 




    TOF = tout(end);
    m_fuel = p.m_dot*TOF;
    dV = (g*p.Isp)*log(p.mass_init/(p.mass_init-m_fuel));
    flip_t = p.t_flip;

    
    h_hat = h/norm(h);

    if dot(h_hat,[0 0 1])>0
        leadvtrail = 1;
    else
        leadvtrail = -1;
    end



end