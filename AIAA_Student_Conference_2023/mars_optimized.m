%mars_optimized


%--Constants for any run--%
MuS=1.32747e20; % Mu for the Sun

rs=149597870700; % m distance from sun to earth

r0=(rs)*[1 0 0].';

v0=sqrt(MuS/norm(r0))*[0 1 0].';  % formula for velocity of a circular 

tperiod=2*pi*sqrt((r0(1,:)^3)/MuS); % formula for one orbit period

d_o= 2.28e+11; %specifes the target orbit

vmt= 24000; %mars orbital speed in m/s

e_m=0;


mass_init=250000%15000; %mass of spacecraft and fuel

%--3g accel---
m_dot=0.053;
thrust=3*mass_init*9.81;
opt_cooef = [51615.043463676346 -0.10108641569407725 0.12944941094690648 ...
             0.56368699337292383 0.028579201156414144 0.014875566842945196 ...
             0.65359358493784514 -0.003066205125677251];

%--1g accel--
m_dot_1=0.028;
thrust_1=mass_init*9.81;
opt_cooef_1=[89403.999090953876 0.20599494657613396 -0.11557025866685769 ...
             -0.15962209667753074 0.70915464496072111 4.6041487744336607 ...
             0.47814755491954941 -0.00053773762753318779];


%--0.5g accel--
m_dot_half=0.02;
thrust_half=0.5*mass_init*9.81;
opt_cooef_half= [126447.07628224068 -0.67883995080904058 2.2641912090027496 ...
             -0.83011366760054539 0.074699585200471652 0.67512669880118481 ...
             0.5881772316169217 -0.0032890612332791766];


mass_init=250000;
thrust=0.019*mass_init*9.81;
m_dot=0.04;



opt_cooef=[664430.999343951	1.04587235276916	0.739785490806324	0.429393100290926	-0.224949960872908	0.100516528969655	1.18976131016172	-0.209860826653552];



