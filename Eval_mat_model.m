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


g=9.81;
MuE = 3.9869044e14;                 % Graviational Parameter for Earth in useful units 
MuM = 4.9048695e12;                 % Gravitational parameter for the Moon
MuS = 1.32712440018e20;              % Graviational parameter for the Sun
mz = 1000.*state(3,:).';

v_em=[moon_dat.moon_vx moon_dat.moon_vy moon_dat.moon_vz].';
[semi_maj,~,~,~,~,~] = posvel2orbitalelements([moon_dat.moon_x(1) moon_dat.moon_y(1) mz(1)].',v_em(:,1),MuM+MuE);
mean_motion = sqrt((MuM+MuE)/semi_maj^3);


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
thrust=366;%mass_init*a0;              % compute the constant thrust given by the initial mass/accel of the spacecraft
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

p.t_flip=0.54*tstop;
th=2.4;

r0=(Re+h)*[cos(th) sin(th) 0].';
v0=(sqrt(MuE/(Re+h)))*[-sin(th) cos(th) 0].';

% % Create ICs for a Lunar Orbit 
% r0 = rf*[1 0 0].' + r_em(:,1);
% v0 = vf*[0 1 0].' + v_em(:,1);

X0 = [r0;v0];
t_span = [0 tstop];

tol = 1e-11;

options = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,X) turn_and_burn_events(t,X,funs,p));

[tout,Xout] = ode113(@(t,X) turn_and_burn_dynamics(t,X,p,funs),t_span,X0,options);


%% Data Processing

% Unpack the geo-centric position and velocites
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

e_vec=(cross(v_s_m(end,:),h))/MuM - r_s_m(end,:)/norm(r_s_m(end,:));

e_s=norm(e_vec); %e of the spacecraft 
%%

figure(1)
hold on
%xlim([0.9*norm(0.5*r_em(:,1)) 1.3*norm(r_em(:,1))])
%ylim([-norm(r_e_m(:,1)) norm(r_e_m(:,1))])
axis equal

%plot(x_data,y_data,"DisplayName","Trajectory" + num2str(k),LineWidth=2)
plot(r_s_e(:,1),r_s_e(:,2),"DisplayName","Trajectory",LineWidth=2)
plot(0,0,'.','MarkerSize',15,'DisplayName','Earth')
plot(r_m_e(:,1),r_m_e(:,2),'DisplayName','Moon')
title("NTP Lunar Trajectory - Tur" + ...
    "n and Burn","FontSize",12)
legend
e_orbit = nsidedpoly(10000, 'Center', [r_m_e(end,1),r_m_e(end,2)], 'Radius',rf);
plot(e_orbit,'FaceColor', 'none', 'EdgeColor','b','LineWidth',2)
% 
%  for i = 1:length(time)
% 
%     h1= plot(x_data(i),y_data(i),'r.','MarkerSize',20);
%     h2= plot(r_emx(i),r_emy(i),'g.','MarkerSize',20);
%     pause(1e-20)
%     if i==length(tout)
%         break
%     end
% 
%     delete(h1)
%     delete(h2)
% end


% figure(2)
% hold on
% plot(Force_earth.time,Force_earth.data,'DisplayName',"Force due to Earth" + num2str(k),LineWidth=3)
% plot(Force_propulsion.time,Force_propulsion.data,'DisplayName',"Force due to Propulsion"+ num2str(k),LineWidth=3)
% plot(Force_moon.time,Force_moon.data,'DisplayName',"Force due to Moon"+ num2str(k),LineWidth=3)
% legend



%% Plot the Trajectories in the rotating frame

x_data = Xout(:,1);
y_data = Xout(:,2);

for i=1:length(tout)
    nu = mean_motion*tout(i);
    r_rot(:,i) = [cos(nu) sin(nu); -sin(nu) cos(nu)]*[x_data(i) y_data(i)].';
    r_moon(:,i) = [cos(nu) sin(nu); -sin(nu) cos(nu)]*[r_m_e(i) r_m_e(i)].';
end



% figure(3)
% plot(r_rot(1,:),r_rot(2,:))
% hold on
% plot(r_moon(1,:),r_moon(2,:),'LineStyle','--')
% plot(r_moon(1,end),r_moon(2,end),'.','MarkerSize',15)



% 
% 
% 
% 
% 
% 
% min_rem = min(r_em_norm);
% if min_rem<Rm
%     error("Collision with Moon")
% end

