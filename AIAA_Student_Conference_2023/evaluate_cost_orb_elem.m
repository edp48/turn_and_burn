%test the cost functon

mars_optimized %if testing new cooef...comment out opt_cooef line in file

t_flip=opt_cooef(1);

a_set=opt_cooef(2:end);

tstop=tperiod; %point at which the model stops running


maxstep=3; %timestep reguator in the simulink model


%set the input to and calls the simulation 
%function to set the coeef in simulink model


sim("orbitEOM_o_e_graph.slx");

%%

%find Simulation end time
time_data=pos.Time;
sizet=size(time_data); %finds the dimentions of time array
time_index=sizet(1);
tf=time_data(time_index); %identifies the last index of the time array (end time)
t_min=tf/60;
t_hrs=t_min/60;
t_days=t_hrs/24;
t_out=[tf t_min t_hrs t_days];

%Final position vector
pos_data=pos.Data;
r=pos_data(end,:);

%Final Velocity Vector
v=vel(end,:);

%Calc ecentricity

h=cross(r,v);

e_vec=(cross(v,h))/MuS - r/norm(r);

e_s=norm(e_vec);


% 
% 
% %percent flip
% 
% 
% disp("radial velocity is "+vr_final)
% disp("tan velocity is "+vtan_final+ " desired is "+ vmt)
% disp("pos is "+r_end+ " desired is "+d_o )
disp("time is (hrs) "+t_hrs )
disp("t_flip is (hrs) "+(t_flip/60)/60)
disp("mass is "+(mass_init-m_dot*tf))
disp("eccentricity is "+e_s)
disp("max speed is "+v_max)

posplot