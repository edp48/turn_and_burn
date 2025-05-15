%position plot script...Universal to both versions of the Matlab Simulink
%Files

figure(1)


plot(pos.Data(:,1),pos.Data(:,2),"Color",'k','LineWidth',2,'DisplayName','Optimal Trajectory')

hold on

e_orbit = nsidedpoly(10000, 'Center', [0,0], 'Radius',149597870700);
plot(e_orbit,'FaceColor', 'none', 'EdgeColor','b','LineWidth',2,'DisplayName','Earth Orbit')

m_orbit= nsidedpoly(10000, 'Center', [0,0], 'Radius',2.28e+11);

plot(m_orbit,'FaceColor', 'none','EdgeColor','r','LineWidth',2,'DisplayName','Mars Orbit')
title('Path of Turn and Burn Trajectory between Mars and Earth Orbits','FontSize',16)
legend('FontSize',13)




figure(2)

plot(pos_norm.time,pos_norm.data,'LineWidth',2)
title('Position vs Time','FontSize',16)
hold on
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Position, m','FontSize',14,'FontWeight','bold')

figure(3)
plot(vel_norm.time,vel_norm.data,'LineWidth',2)
title('Velocity vs Time','FontSize',16)
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Velocity, m/s','FontSize',14,'FontWeight','bold')

figure(4)
plot(accel_norm.time,accel_norm.data,'LineWidth',2)
title('Acceleration vs Time','FontSize',16)
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Acceleration, m/s^2','FontSize',14,'FontWeight','bold')

figure(5)
plot(pos_norm.time,t_angle,'LineWidth',2)
title('Angle of Thrust relative to Major Axis','FontSize',16)
xlabel('Time, s','FontSize',14,'FontWeight','bold')
ylabel('Angle of T Vector and Horizontal Axis, rad','FontSize',14,'FontWeight','bold')
