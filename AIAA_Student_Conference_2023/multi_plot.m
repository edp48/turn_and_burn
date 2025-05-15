%% Plot multiple Acceleration Conditions


close all
clear
clc

mars_optimized %if testing new cooef...comment out opt_cooef line in file

tstop=tperiod; %point at which the model stops running

maxstep=1;








%% Graph 1/2 g

t_flip=opt_cooef_half(1);

a_set=opt_cooef_half(2:end);
m_dot=m_dot_half;
thrust=thrust_half;

sim("orbitEOM_o_e_graph.slx");

figure(1)
plot(pos.Data(:,1),pos.Data(:,2),'LineWidth',2,'DisplayName','g/2')
title('Path of Turn and Burn Trajectory between Mars and Earth Orbits','FontSize',16)
legend('FontSize',13)
hold on



figure(2)
plot(pos_norm.time,pos_norm.data,'LineWidth',2,'DisplayName','g/2')
title('Position vs Time','FontSize',16)
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Position, m','FontSize',14,'FontWeight','bold')
hold on
legend('FontSize',13)


figure(3)
plot(vel_norm.time,vel_norm.data,'LineWidth',2,'DisplayName','g/2')
title('Velocity vs Time','FontSize',16)
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Velocity, m/s','FontSize',14,'FontWeight','bold')
hold on
legend('FontSize',13)

figure(4)
plot(accel_norm.time,accel_norm.data,'LineWidth',2,'DisplayName','g/2')
title('Acceleration vs Time','FontSize',16)
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Acceleration, m/s^2','FontSize',14,'FontWeight','bold')
hold on
legend('FontSize',13)

figure(5)
plot(pos_norm.time,t_angle,'LineWidth',2,'DisplayName','g/2')
title('Angle of Thrust relative to Major Axis','FontSize',16)
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Angle of T Vector and Horizontal Axis, rad','FontSize',14,'FontWeight','bold')
hold on
legend('FontSize',13)

%% Graph 1g

t_flip=opt_cooef_1(1);

a_set=opt_cooef_1(2:end);
m_dot=m_dot_1;
thrust=thrust_1;

sim("orbitEOM_o_e_graph.slx");

figure(1)
plot(pos.Data(:,1),pos.Data(:,2),"Color",'k','LineWidth',2,'DisplayName','1g')


figure(2)
plot(pos_norm.time,pos_norm.data,'LineWidth',2,'DisplayName','1g')

figure(3)
plot(vel_norm.time,vel_norm.data,'LineWidth',2,'DisplayName','1g')


figure(4)
plot(accel_norm.time,accel_norm.data,'LineWidth',2,'DisplayName','1g')


figure(5)
plot(pos_norm.time,t_angle,'LineWidth',2,'DisplayName','1g')


%% 3g Plots

t_flip=opt_cooef_3(1);

a_set=opt_cooef_3(2:end);
m_dot=m_dot_3;
thrust=thrust_3;

sim("orbitEOM_o_e_graph.slx");

figure(1)
plot(pos.Data(:,1),pos.Data(:,2),'LineWidth',2,'DisplayName','3g')
e_orbit = nsidedpoly(10000, 'Center', [0,0], 'Radius',149597870700);
plot(e_orbit,'FaceColor', 'none', 'EdgeColor','b','LineWidth',2)
m_orbit= nsidedpoly(10000, 'Center', [0,0], 'Radius',2.28e+11);
plot(m_orbit,'FaceColor', 'none','EdgeColor','r','LineWidth',2)



figure(2)
plot(pos_norm.time,pos_norm.data,'LineWidth',2,'DisplayName','3g')

figure(3)
plot(vel_norm.time,vel_norm.data,'LineWidth',2,'DisplayName','3g')

figure(4)
plot(accel_norm.time,accel_norm.data,'LineWidth',2,'DisplayName','3g')

figure(5)
plot(pos_norm.time,t_angle,'LineWidth',2,'DisplayName','3g')




