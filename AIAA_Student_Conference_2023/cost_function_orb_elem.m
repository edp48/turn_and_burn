%cost_function 

function j= cost_function_orb_elem(guess_init)

    MuS=1.32747e20;
    d_o= 2.28e+11;
    e_m=0;

    %set the input to and calls the simulation 

    t_set=guess_init(1);

    cooef_set=guess_init(2:end);
    
    assignin("base","a_set",cooef_set); %function to set the coeef in simulink model
    
    assignin("base","t_flip",t_set);

    sim("orbitEOM_orbital_elements.slx");
    
    %Final velocity vector
    pos_data=pos.Data;
    r=pos_data(end,:);
    
    %Final Velocity Vector
    v=vel(end,:);
    
    %Calc ecentricity
    
    h=cross(r,v);
    
    e_vec=(cross(v,h))/MuS - r/norm(r);
    
    e_s=norm(e_vec); %e of the spacecraft 
    
    c=[0.000000000001 1 1e30];
    
    tf=pos.Time;
    
    %Cost Function
    j=c(1)*(tf).^2 + c(2)*(norm(r)-d_o).^2 + c(3)*(e_s-e_m)^2;
end