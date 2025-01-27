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

    dXdt(4:6) = (1/m_time)*((-m_time)*((p.MuE/(norm(r_s_e)^3))*r_s_e + (p.MuM/(norm(r_s_m)^3))*r_s_m) + thrust*v_hat);


   

   end