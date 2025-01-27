function [value,isterminal,direction] = turn_and_burn_events(t,X,funs,p)

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




    v_m_e(1) = funs.Intvx(t);
    v_m_e(2) = funs.Intvy(t);
    v_m_e(3) = funs.Intvz(t);
    
    % Calculate the spacecraft's lunar-centric pos/vel
    %r_s_m = r_s_e-r_m_e;
    v_s_m = v_s_e-v_m_e;
    
    
    %Calculate eccentricity about the moon
    h=cross(r_s_m(end,:),v_s_m(end,:));
    
    e_vec2=(cross(v_s_m(end,:),h))/p.MuM - r_s_m(end,:)/norm(r_s_m(end,:));
    
    e_s2=norm(e_vec2); %e of the spacecraft 


value = e_s2-0.6;
isterminal = 1;
direction = 0;




end