% Used to Recreate figures and Data presented in "Optimal Cislunar
% Trajectories with Continuous, High-Thrust Nuclear-Thermal Propulsion" at
% Scitech 2025. Consult README before running

clear
clc
close all

% Use to determine if figs are saved as files in current directory
% Write as true or false
savefig = false;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%      Load In the .MAT file and create the other relevant variables 
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


load Exhaust_prop_data.mat

% Create the Moon's State with MICE
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

r_em=[moon_dat.moon_x moon_dat.moon_y moon_dat.moon_z].'; 


% Generate Interpolation Functions for use later 
funs.Intx = griddedInterpolant(moon_dat.time,moon_dat.moon_x);
funs.Inty = griddedInterpolant(moon_dat.time,moon_dat.moon_y);
funs.Intz = griddedInterpolant(moon_dat.time,moon_dat.moon_z);
funs.Intvx = griddedInterpolant(moon_dat.time,moon_dat.moon_vx);
funs.Intvy = griddedInterpolant(moon_dat.time,moon_dat.moon_vy);
funs.Intvz = griddedInterpolant(moon_dat.time,moon_dat.moon_vz);

%Generate Fixed Parameters {Known constants, Spacecraft and Orbital}
g=9.81;
MuE = 3.9869044e14;                
MuM = 4.9048695e12;                 
MuS = 1.32712440018e20;             

%Spacecraft Parameters
mass_init=20000;
mf=9000; 
isp=1000;     


Re=6378176;                         % Radius of the earth [m]
Rm=1740000;                         % Radius of the moon 
h=35000000;                         % GEO orbit [m]
h_lunar = 100000;                   % LLO orbit [m]


% Package Parameters for use later
p.Rm = Rm;
p.h_lunar = h_lunar;
p.Re = Re;
p.h = h;
p.MuE = MuE;
p.MuM = MuM;
p.mass_init = mass_init;
p.Isp = isp;

% Pre-create Vectors for Speed
TOF_v = zeros(size(thrust_v));
burned_prop = zeros(size(thrust_v));
dV_v = zeros(size(thrust_v));
flip_t_v = zeros(size(thrust_v));
eccent = zeros(size(thrust_v));


% Thrust Dependent Parameters
for t_ind = 1:length(thrust_v)

    thrust=thrust_v(t_ind);             
    m_dot=thrust/(g*isp);
    tstop=mf/m_dot;
    
 
    p.thrust = thrust;
    p.m_dot = m_dot;
    p.tstop = tstop;

    % Extract the relevent quantites from the Simulation
    [TOF_v(t_ind),burned_prop(t_ind),dV_v(t_ind),flip_t_v(t_ind),...
        eccent(t_ind),~] = extract_info(opt_initial_vec(:,t_ind),p,funs);

end

% Some final Data Manipulation
TOF_days = TOF_v/86400;
flip_perc = flip_t_v./TOF_v;




%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                          Generate the Needed Figures
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%Some Plotting Parameters for uniformity
ftsz = 18;
mkrsz = 50;



%%%%%%%%%%%%%%%%%%%%%%%%%%% Time of Flight Figure %%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
box on

% Data
scatter(thrust_v,TOF_days,mkrsz,[0.6350 0.0780 0.1840],"^","filled",...
    'DisplayName',"Leading Flybys",'LineWidth',1)
plot(thrust_v,3*ones(size(thrust_v)),"LineStyle","--","LineWidth",2,...
    "DisplayName","Lunar Free Return")

% Labels
xlabel("Thrust [N]",'Interpreter','latex','FontSize',ftsz)
ylabel("Time of Flight [days]", 'Interpreter','latex',FontSize=ftsz)
xlim([290 410])

% Legand and Axis
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnd = legend("Interpreter","latex",FontSize=ftsz);
lgnd.Location = "southwest";



if savefig
    print(gcf, 'min_TOF', '-depsc', '-r300');
    print(gcf, 'min_TOF', '-dsvg');
end


%%%%%%%%%%%%%%%%%%% Propellent Mass Burned Figure %%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
box on

% Data
scatter(thrust_v,burned_prop,mkrsz,[0.6350 0.0780 0.1840],"^","filled",...
    'DisplayName',"Leading Flybys",'LineWidth',1)

% Labels
xlabel("Thrust [N]",'Interpreter','latex','FontSize',ftsz)
ylabel("Propellant Mass Expelled [kg]", 'Interpreter','latex',FontSize=ftsz)
xlim([290 410])

% Legand and Axis
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnd = legend("Interpreter","latex",FontSize=ftsz);
lgnd.Location = "southeast";

if savefig
    print(gcf, 'min_Prop_Mass', '-depsc', '-r300');
    print(gcf, 'min_Prop_Mass', '-dsvg');
end 


%%%%%%%%%%%%%%%%%%%%%%%%%% Delta V Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
box on

% Data
scatter(thrust_v,dV_v./1000,mkrsz,[0.6350 0.0780 0.1840],"^","filled",...
    'DisplayName',"Leading Flybys",'LineWidth',1)
plot(thrust_v,3*ones(size(thrust_v)),"LineStyle","--","LineWidth",2,...
    "DisplayName","Lunar Free Return")

% Labels
xlabel("Thrust [N]",'Interpreter','latex','FontSize',ftsz)
ylabel("Maneuver $\Delta V$, [$\frac{km}{s}$]", 'Interpreter',...
    'latex',FontSize=ftsz)
xlim([290 410])
ylim([2.9 6])

% Legand and Axis
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnddv = legend("Interpreter","latex",FontSize=ftsz);
lgnddv.Position = [0.5520    0.7448    0.3373    0.1136];

if savefig
    print(gcf, 'min_dV', '-depsc', '-r300');
    print(gcf, 'min_dV', '-dsvg');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Flip Time Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the trendline
ff = polyfit(thrust_v,flip_perc,1);
bestfit = ff(1)*thrust_v + ff(2);


figure(4)
hold on
box on

% Data
plot(thrust_v,bestfit,'LineWidth',4,"DisplayName",...
    ['3.7396E-04 $\times T$ + ',num2str(ff(2))],"Color",'#3FD5EA')
scatter(thrust_v,flip_perc,mkrsz,[0.6350 0.0780 0.1840],"^","filled",...
    'DisplayName',"Leading Flybys",'LineWidth',1)


% Labels
xlabel("Thrust [N]",'Interpreter','latex','FontSize',ftsz)
ylabel("Flip Time [\% of $\tau$]", 'Interpreter','latex',FontSize=ftsz)
xlim([290 410])

% Legand and Axis
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnd = legend("Interpreter","latex",FontSize=ftsz);
lgnd.Location = "southeast";

if savefig
    print(gcf, 'min_flip_time', '-depsc', '-r300');
    print(gcf, 'min_flip_time', '-dsvg');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Eccentricity Figure %%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
hold on
box on


% Data
scatter(thrust_v,eccent,mkrsz,[0.6350 0.0780 0.1840],"^",...
    "filled",'DisplayName',"Leading Flybys",'LineWidth',1)

% Labels
xlabel("Thrust [N]",'Interpreter','latex','FontSize',ftsz)
ylabel("Eccentricity at $\tau$",'Interpreter','latex',FontSize=ftsz)
ylim([0 1])
xlim([290 410])

% Legand and Axis
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
legend("Interpreter","latex",FontSize=ftsz);


if savefig
    print(gcf, 'min_eccent', '-depsc', '-r300');
    print(gcf, 'min_eccent', '-dsvg');
end


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                     Generate the Bounding Trajectories
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


%%%%%%%%%%%%%%%%%%%%%%%%%% Lowest Thrust Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


thrust = thrust_v(1);
m_dot=thrust/(g*isp);
tstop=mf/m_dot;


p.thrust = thrust;
p.m_dot = m_dot;
p.tstop = tstop;

p.t_flip=opt_initial_vec(1,1);
th=opt_initial_vec(2,1);

r0=(Re+h)*[cos(th) sin(th) 0].';
v0=(sqrt(MuE/(Re+h)))*[-sin(th) cos(th) 0].';


X0 = [r0;v0];
t_span = [0 tstop];


tol = 1e-11;

options = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,X) turn_and_burn_events(t,X,funs,p));

[toutmin,Xoutmin] = ode113(@(t,X) turn_and_burn_dynamics(t,X,p,funs),t_span,X0,options);


%%%%%%%%%%%%%%%%%%%%%%%%% Largest Thrust Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


thrust=thrust_v(end);
m_dot=thrust/(g*isp);
tstop=mf/m_dot;

p.thrust = thrust;
p.m_dot = m_dot;
p.tstop = tstop;

 
p.t_flip=opt_initial_vec(1,end);
th=opt_initial_vec(2,end);

r0=(Re+h)*[cos(th) sin(th) 0].';
v0=(sqrt(MuE/(Re+h)))*[-sin(th) cos(th) 0].';


X0 = [r0;v0];
t_span = 0:5:tstop;


tol = 3e-14;

options = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,X) turn_and_burn_events(t,X,funs,p));

[toutmax,Xoutmax] = ode113(@(t,X) turn_and_burn_dynamics(t,X,p,funs),t_span,X0,options);




%%%%%%%%%%%%%%%%%%%%%%%%%% Data Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Min Data Extraction

DU_con = 384400000;

r_s_e_min = Xoutmin(:,1:3)./DU_con;


r_m_emin(:,1) = funs.Intx(toutmin)./DU_con;
r_m_emin(:,2) = funs.Inty(toutmin)./DU_con;
r_m_emin(:,3) = funs.Intz(toutmin)./DU_con;

%Max Data Extraction

r_s_e_max = Xoutmax(:,1:3)./DU_con;


r_m_emax(:,1) = funs.Intx(toutmax)./DU_con;
r_m_emax(:,2) = funs.Inty(toutmax)./DU_con;
r_m_emax(:,3) = funs.Intz(toutmax)./DU_con;


% Create the plot
figure(6)
hold on
axis equal
ylim([-0.2 0.72])
plot(r_s_e_min(:,1),r_s_e_min(:,2),"DisplayName","Trajectory..." + ...
    " Lower Bound",LineWidth=2)

plot(r_m_emin(end,1),r_m_emin(end,2),'.','MarkerSize',25,'DisplayName'...
    ,'Moon Position Lower Bound')

plot(r_s_e_max(:,1),r_s_e_max(:,2),"DisplayName",...
    "Trajectory Upper Bound",LineWidth=2)
plot(0,0,'.','MarkerSize',30,'DisplayName','Earth')
plot(r_m_emax(end,1),r_m_emax(end,2),'.','MarkerSize',25,...
    'DisplayName','Moon Position Upper Bound')
xlabel("X-position [DU]","Interpreter","latex",FontSize=ftsz)
ylabel("Y-position [DU]","Interpreter","latex",FontSize=ftsz)


ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnd = legend("Interpreter","latex",FontSize=ftsz);
lgnd.Location = "northwest";
box on


if savefig
    print(gcf, 'min_Bounded_Trajectories', '-depsc', '-r300');
    print(gcf, 'min_Bounded_Trajectories', '-dsvg');
end 






%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                     Functions used in this File
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




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

    dXdt(4:6) = (1/m_time)*((-m_time)*((p.MuE/(norm(r_s_e)^3))*r_s_e + ...
        (p.MuM/(norm(r_s_m)^3))*r_s_m) + thrust*v_hat);


   

end

function [TOF,m_fuel,dV,flip_t,e_s,leadvtrail] = extract_info(opt_inputs,p,funs)

    % Run the Simulation
    g=9.81;
    MuE = 3.9869044e14; 
    %Orbital Parameters
    Re=6378176;                         % Radius of the earth [m]
    %h=400000; %ALSO CHANGE IN COST
    h=35000000;% LEO orbit [m]
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
