function [r,v] = orbitalelements2posvel(a,e,Omega,I,omega,mu,nu)


r_p_mag = (a*(1-e^2))/(1+e*cos(nu));

r_p = [r_p_mag*cos(nu) r_p_mag*sin(nu) 0].';


l=a*(1-e^2);

v_p = sqrt(mu/l)*[-sin(nu) e+cos(nu) 0].';

rot1=[cos(Omega) sin(Omega) 0;-sin(Omega) cos(Omega) 0; 0 0 1];

rot2=[1 0 0; 0 cos(I) sin(I); 0 -sin(I) cos(I)];

rot3=[cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

PCI=rot3*rot2*rot1;

ICP=PCI.';

r=ICP*r_p;
v=ICP*v_p;

end