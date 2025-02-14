% Used to Recreate figures and Data presented in "Optimal Cislunar
% Trajectories with Continuous, High-Thrust Nuclear-Thermal Propulsion" at
% Scitech 2025. Consult README before running

clear
clc
close all

savefig = false;

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
% Sim Data
    
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

p.t_flip=148205.570730942;
th = 2.01254733806579;


r0=(Re+h)*[cos(th) sin(th) 0].';
v0=(sqrt(MuE/(Re+h)))*[-sin(th) cos(th) 0].';


X0 = [r0;v0];
t_span = [0 tstop];

tol = 1e-11;

options = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,X) turn_and_burn_events(t,X,funs,p));

[tout,Xout] = ode113(@(t,X) turn_and_burn_dynamics(t,X,p,funs),t_span,X0,options);


%% Data Processing


[accel_thrust,accel_moon,accel_earth]=plot_accel(Xout,tout,p,funs);


ftsz=18;

figure(1)
hold on
box on

plot(tout,abs(accel_thrust),'LineWidth',2,'DisplayName',"Thrust")
plot(tout,accel_earth,'LineWidth',2,'DisplayName',"Earth")
plot(tout,accel_moon,'LineWidth',2,'DisplayName',"Moon")
xline(p.t_flip,'--','LineWidth',2,"DisplayName","$\tau_m$")
xlabel("Time, s","FontSize",ftsz,"Interpreter","latex")
ylabel("Acceleration, $\frac{m}{s^2}$","FontSize",ftsz,"Interpreter","latex")

ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnd = legend("Interpreter","latex",FontSize=ftsz);
lgnd.Location = "northeast";


if savefig
    print(gcf, 'accel_components', '-depsc', '-r300');
    print(gcf, 'accel_components', '-dsvg');
end







%% Functions

function [accel_thrust,accel_moon,accel_earth]=plot_accel(Xout,tout,p,funs)

% Pre-load variables for speed
accel_thrust = zeros(size(tout));
accel_moon = zeros(size(tout));
accel_earth = zeros(size(tout));


for a=1:length(tout)

    X=Xout(a,:);
    t=tout(a);

    r_s_e = X(1:3);
    r_s_e = reshape(r_s_e,1,3);

   
    
    r_m_e(1) = funs.Intx(t);
    r_m_e(2) = funs.Inty(t);
    r_m_e(3) = funs.Intz(t);
    r_m_e = reshape(r_m_e,1,3);
    

    r_s_m = r_s_e - r_m_e;


    
    m_time = p.mass_init - p.m_dot*t;
    
    accel_thrust(a) = (1/m_time)*p.thrust*((p.t_flip-t)/norm(p.t_flip-t));
    accel_moon(a) = norm((p.MuM/(norm(r_s_m)^3))*r_s_m);
    accel_earth(a) = norm((p.MuE/(norm(r_s_e)^3))*r_s_e);

end






end





function [value,isterminal,direction] = turn_and_burn_events(t,X,funs,p)

r_s_e = X(1:3);
    r_s_e = reshape(r_s_e,1,3);
    v_s_e = X(4:6);
    v_s_e = reshape(v_s_e,1,3);
   
    

    
    r_m_e(1) = funs.Intx(t);
    r_m_e(2) = funs.Inty(t);
    r_m_e(3) = funs.Intz(t);
    r_m_e = reshape(r_m_e,1,3);
    

    r_s_m = r_s_e - r_m_e;




    v_m_e(1) = funs.Intvx(t);
    v_m_e(2) = funs.Intvy(t);
    v_m_e(3) = funs.Intvz(t);
    
    % Calculate the spacecraft's lunar-centric pos/vel
    %r_s_m = r_s_e-r_m_e;
    v_s_m = v_s_e-v_m_e;
    
    
    %Calculate eccentricity about the moon
    h=cross(r_s_m(end,:),v_s_m(end,:));
    
    e_vec2=(cross(v_s_m(end,:),h))/p.MuM - r_s_m(end,:)/norm(r_s_m(end,:));
    
    e_s2=norm(e_vec2); %e of the spacecraft 


value = e_s2-0.0000006;
isterminal = 1;
direction = 0;




end

function [dXdt] = turn_and_burn_dynamics(t,X,p,funs)
    
    r_s_e = X(1:3);
    r_s_e = reshape(r_s_e,1,3);
    v_s_e = X(4:6);
    v_s_e = reshape(v_s_e,1,3);
    v_hat = v_s_e/(norm(v_s_e) +0.2);
    

    
    r_m_e(1) = funs.Intx(t);
    r_m_e(2) = funs.Inty(t);
    r_m_e(3) = funs.Intz(t);
    r_m_e = reshape(r_m_e,1,3);
    

    r_s_m = r_s_e - r_m_e;


    dXdt = zeros(6,1);

    dXdt(1) = X(4);
    dXdt(2) = X(5);
    dXdt(3) = X(6);
    
    m_time = p.mass_init - p.m_dot*t;
    
    thrust = p.thrust*((p.t_flip-t)/norm(p.t_flip-t));

    dXdt(4:6) = (1/m_time)*((-m_time)*((p.MuE/(norm(r_s_e)^3))*r_s_e + (p.MuM/(norm(r_s_m)^3))*r_s_m) + thrust*v_hat);


   

   end
