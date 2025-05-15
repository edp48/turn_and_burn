function [dXdt] = Dynamic_Model_Lambert(X,p)

   
    
    r = X(1:3,:);
    v = X(4:6,:);


    dXdt = zeros(size(X));

    dXdt(1,:) = X(4,:);
    dXdt(2,:) = X(5,:);
    dXdt(3,:) = X(6,:);
    

    dXdt(4:6,:) = -(p.mu./(vecnorm(r).^3)).*r;
end