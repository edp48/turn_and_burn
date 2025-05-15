function [dXdt] = Dynamic_Model(t,X,ang,p)

   
    
    r = X(1:3,:);
    v = X(4:6,:);
    
    
    v_hat = v./vecnorm(v);

    t_hat = zeros(size(v_hat));

    for k = 1:length(t)
        t_hat(:,k) = rot_z(v_hat(:,k),ang(k));
    end
    

    dXdt = zeros(size(X));

    dXdt(1,:) = X(4,:);
    dXdt(2,:) = X(5,:);
    dXdt(3,:) = X(6,:);
    
    m_time = p.mass_init - p.m_dot.*t;

    thrust = p.thrust;

    dXdt(4:6,:) = (1./m_time).*((-m_time).*((p.mu./(vecnorm(r).^3)).*r) + thrust.*t_hat);






    function vec_out = rot_z(vec_in,theta_in)
        
        sz = size(vec_in);
        if sz ~= [3,1]
            error("Wrong Size Vector")
        end 
        
        theta = theta_in; %deg2rad(theta_in);
        
        R_z = [cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];
        
        vec_out = R_z*vec_in;
    
    
    end 

   

   end