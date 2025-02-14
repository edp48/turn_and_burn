% Used to Recreate figures and Data presented in "Optimal Cislunar
% Trajectories with Continuous, High-Thrust Nuclear-Thermal Propulsion" at
% Scitech 2025. Consult README before running

clear
close all
clc

savefig = false;

%Paramaters for the model
re=6378000; %m (radius of the earth)
rm=1737400; %m
mue=3.986e14; %Mu for the earth
mum=4.9028e12; %Mu for the moon

mass_sun=1.9891e30; %kg
mass_earth=5.97e24; %kg
mass_moon=7.348e22; %kg
g=9.81;

r_s_e=149597870700; %m
r_e_m=384400000; %m

%Establish the starting and ending orbits of the spacecraft
R_SOI_e=r_s_e*(mass_earth/mass_sun)^(2/5); %Earth and Sun Sphere of Influance
R_SOI_m=r_e_m*(mass_moon/mass_earth)^(2/5); % Moon and Earth Spehere of Influance

% Boundry Conditions: vesc at both SOI
v0=sqrt(2*mue/R_SOI_e);
vf=sqrt(2*mum/R_SOI_m);

%distance of trajectory (from earth to moon)
r_em=384467000; %m


%establish the mass of spacecraft and fuel
md=130000; %kg (dry mass of the spacecraft)
mf=70000; %kg (mass of fuel)
mtot=md+mf;



isp=1000;

%Calcualte the mdot for a given accel
a0=.0014*9.81; %m/s^2 initial accel of the spacecraft...used to find mdot
thrust=(mtot)*a0; %compute the constant thrust given by the initial mass/accel of the spacecraft
mdot=thrust/(g*isp);

% calculate the end time for the manuver, assuming use of all
% propellant
tf=mf/mdot;

%This is the new turn time generalized for some V0 and Vf as defined by
%starting and ending orbits
turn_time=(1/mdot)*(mtot-sqrt((mtot^2 -mtot*mdot*tf)*exp(-(vf-v0)/(g*isp))));

%given that we will use all the alloted fuel, the two time spans are: 
t1=0:0.5:turn_time;
t2=turn_time:0.5:tf;


%pack constants for use in ODE45
p=[];
p.g=g;
p.isp=isp;
p.md=md;
p.mdot=mdot;
p.mtot=mtot;
p.v0=v0;
p.turn_time=turn_time;




%ODE45 for the first half
r0_ODE=0; %Only IC needed for integration is starting position (0 for linear analysis)

f_first =@(t1,zin) Velocitypt1(t1,zin,p);
[t_first,x_first]=ode45(f_first,t1,r0_ODE);

%ODE45 for second half
r2_0_ODE=x_first(end); %IC for pt 2 based off the final position of pt 1


f_last =@(t2,zin) Velocitypt2(t2,zin,p);
[t_last,x_last]=ode45(f_last,t2,r2_0_ODE);


%% Graph everything
ftsz=18;

v1=isp.*g.*log(mtot./(mtot-mdot.*t1)) + v0;
v2=isp.*g.*log((mtot.*(mtot-mdot.*t2))/((mtot-mdot*turn_time)^2)) + v0;

a1=(mdot.*g.*isp)./(mtot-mdot.*t1);
a2=-(mdot.*g.*isp)./(mtot-mdot.*t2);

tnorm = max(t2);
a_norm = max(max(abs(a1)),max(abs(a2)));
v_norm = max(max(v1),max(v2));
p_norm = max(max(x_first),max(x_last));


figure()
hold on
plot(t1/tnorm,x_first/p_norm,'displayname','Before Flip',LineWidth=3)
plot(t2/tnorm,x_last/p_norm,'displayname','After Flip',LineWidth=3)
xlabel("Time $\frac{t}{\tau}$","Interpreter","latex",'FontSize',ftsz,'FontWeight','bold')
ylabel("Normalized Position","Interpreter","latex",'FontSize',ftsz,'FontWeight','bold')
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnd = legend("Interpreter","latex",FontSize=ftsz);
lgnd.Location = "southeast";
box on

if savefig
    print(gcf, 'pos_closed_form', '-depsc', '-r300');
    print(gcf, 'pos_closed_form', '-dsvg');
end


figure()
hold on
box on
plot(t1/tnorm,v1/v_norm,'displayname','Before Flip',LineWidth=3)
plot(t2/tnorm,v2/v_norm,'displayname','After Flip',LineWidth=3)
xlabel("Time $\frac{t}{\tau}$","Interpreter","latex",'FontSize',ftsz,'FontWeight','bold')
ylabel("Normalized Velocity","Interpreter","latex",'FontSize',ftsz,'FontWeight','bold')
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
legend("Interpreter","latex",FontSize=ftsz);

if savefig
    print(gcf, 'vel_closed_form', '-depsc', '-r300');
    print(gcf, 'vel_closed_form', '-dsvg');
end 


figure()
hold on
box on
plot(t1/tnorm,a1/a_norm,'displayname','Before Flip',LineWidth=3)
plot(t2/tnorm,a2/a_norm,'displayname','After Flip',LineWidth=3)
xlabel("Time $\frac{t}{\tau}$","Interpreter","latex",'FontSize',ftsz,'FontWeight','bold')
ylabel("Normalized Acceleration","Interpreter","latex",'FontSize',ftsz,'FontWeight','bold')
ax = gca;
ax.FontSize = ftsz;
ax.TickLabelInterpreter = "latex";
lgnd = legend("Interpreter","latex",FontSize=ftsz);

if savefig
    print(gcf, 'accel_closed_form', '-depsc', '-r300');
    print(gcf, 'accel_closed_form', '-dsvg');
end



%% Functions

%function for the first burn
function dz=Velocitypt1(t1,~,p)
    %Unpack Parameters
    isp=p.isp;
    g=p.g;
    mdot=p.mdot;
    mtot=p.mtot;
    v0=p.v0;

    %Erik's generalization of Daniel's work: Velocity first half
    v1=isp.*g.*log(mtot./(mtot-mdot.*t1)) + v0;


    dz=v1;

end

%function for the second burn
function dz=Velocitypt2(t2,~,p)
    %Unpack parameters
    isp=p.isp;
    g=p.g;
    mdot=p.mdot;
    mtot=p.mtot;
    v0=p.v0;
    turn_time=p.turn_time;

    %Erik's generalization of Daniel's work: Velocity second half
    v2=isp.*g.*log((mtot*(mtot-mdot*t2))/((mtot-mdot*turn_time)^2)) + v0;


    dz=v2;

end