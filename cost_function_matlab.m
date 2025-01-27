%% Cost Function

function j=cost_function_matlab(guess_init,p,funs)

    rf=p.Rm+p.h_lunar;

    %set the input to and calls the simulation 

    p.t_flip=guess_init(1);

    angle_set=guess_init(2);

    r0=(p.Re+p.h)*[cos(angle_set) sin(angle_set) 0].';
    v0=(sqrt(p.MuE/(p.Re+p.h)))*[-sin(angle_set) cos(angle_set) 0].';
        
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


    c=[0.000 1e3 1e4];

    tf=tout(end);

    %Cost Function
    %j=c(1)*(tf).^2 + c(2)*(norm(r)-rf).^2 + c(3)*(e_s-e_m)^2;
    j=c(1)*(tf).^2 + c(2)*10/(1+exp(-(19*norm(r_s_m)-(rf)))).^2 + c(3)*10/(1+exp(-(100*e_s-75)));

%23


end 



