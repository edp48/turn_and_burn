

mars_set_up
numb=1;
iter_numb=3;

%

for n=1:iter_numb


maxstep=3; %timestep reguator in the simulink model
tstop=tperiod/2;



%make changes to the options for fmincon to make it less "weak"
% 
options = optimset('MaxFunEval',10,'MaxIter',10,'Display','final','TolFun',1e-30,'TolX',1e-30);

[opt_cooef,fval,exitflag,output]= fminsearch(@cost_function_orb_elem,guess_init,options);

disp("Iteration #"+ numb+  ": cost value is " + fval)

guess_init=opt_cooef;
numb=numb+1;


end


























% %do_optimization
% 
% %set Global Variables
% %global MuS h m r0 v0 tperiod tstop maxstep d_o vmt thrust 
% 
% mars_set_up
% 
% %
% 
% mass_init=15000; %mass of spacecraft and fuel
% 
% maxstep=5; %timestep reguator in the simulink model
% 
% % %1 g conditions
% % t_0=89409;
% % thrust=*mass_init*9.81;
% % tstop=
% 
% %3 g conditions
% t_0=51700;
% m_dot=0.053; %0.055
% thrust=3*mass_init*9.81;
% tstop=tperiod/2;
% 
% % %10 g conditions
% % t_0=28266;
% %m_dot=0.08;
% %thrust=10*mass_init*9.81;
% %tstop=187000
% %guess_init=[t_0 0 0 0 0 0 0 0]; %initial guess of cooef.
% 
% % % rng default % For reproducibility
% % % gs = GlobalSearch('Display','iter');
% % problem = createOptimProblem('fmincon','x0',guess_init,'objective',@cost_function_flip);
% % ms=MultiStart('Display','final');
% % opt_cooef = run(ms,problem,500);
% 
% 
% %options = optimoptions('fmincon',StepTolerance=1e-20, OptimalityTolerance=1e-20);
% % 
% % 
% % A = [];
% % b = [];
% % Aeq = [];
% % beq = [];
% % lb = [];
% % ub = [];
% % nonlcon = [];
% % 
% % 
% %[opt_cooef,fval,exitflag,output]= fmincon(@cost_function_flip,guess_init,A,b,Aeq,beq,lb,ub,nonlcon,options);
% 
% %make changes to the options for fmincon to make it less "weak"
% % 
% options = optimset('MaxFunEval',7000,'MaxIter',5000,'Display','iter','TolFun',1e-300,'TolX',1e-300);
% 
% [opt_cooef,fval,exitflag,output]= fminsearch(@cost_function_orb_elem,guess_init,options);



