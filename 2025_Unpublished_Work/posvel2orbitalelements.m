function [a,e,Omega,I,omega,nu] = posvel2orbitalelements(r,v,mu)

%Establish unit vectors
e1 = [1 0 0].';
e2 = [0 1 0].';
e3 = [0 0 1].';

r_hat = r/norm(r);
v_hat = v/norm(v);

%Right ascension of the assending node
h_v = cross(r,v);

h_hat = h_v/norm(h_v);


n_hat = cross(e3,h_v)/norm(cross(e3,h_v)); %unit vector line of nodes

%handles the special casw where I=0
if abs(dot(h_hat,e3))==1
    n_hat = e1;
    %disp('equitorial')
end

if dot(n_hat,e2) > 0 % to determine if Omega is (0,pi) or (pi,2pi)
    Omega = acos(dot(n_hat,e1));
else
    Omega = 2*pi-acos(dot(n_hat,e1));
end

% Eccentricity
e_v = cross(v,h_v)./mu - r/norm(r);

e = norm(e_v);

e_tol = 1e-12;

e_hat = e_v/e;

if e<e_tol
    e_hat = r_hat;
    %disp('circular')
end

%True Anomaly
if dot(r,v)>0
    t = dot(e_hat,r_hat);
    if abs(t)>1 && (abs(t)-1)<1e-14
        t=1;
        
    end
    nu= acos(t);
else 
    t = dot(e_hat,r_hat);
    if abs(t)>1 && (abs(t)-1)<1e-14
        t=1;
    end
    nu = 2*pi - acos(t);
end

%Inclination
I=acos(dot(h_hat,e3));


%argument of periapsis
d=cross(h_hat,n_hat);

if dot(e_hat,e3) >=0
    omega = acos(dot(e_hat,n_hat));
else
    omega = 2*pi-acos(dot(e_hat,n_hat));
end

%Circular Orbit Case

if e < (1-e_tol)
    a=norm(r)*((1+ e*cos(nu))/(1-e^2));
    %disp("Closed Orbit")
end


%Parabolic Case

if (e >(1-e_tol))&&(e<(1+e_tol))
    l=norm(r)*(1+e*cos(nu));
    a=l/2;
    %disp("Parabolic Orbit")
end

%Hyperbolic Case

if e>1+e_tol
    tan_nu_2 = tan(nu/2);
    e_express = sqrt((e+1)/(e-1));
    H=2*atanh(tan_nu_2/e_express);
    a=norm(r)/(1-e*cosh(H));
    %disp('Hyperbolic Orbit')
end