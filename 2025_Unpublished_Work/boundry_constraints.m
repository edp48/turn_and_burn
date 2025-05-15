function [C_ineq,C_eq] = boundry_constraints(X0,Xf,p)



%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Initial BCs $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

r0 = X0(1:3);
v0 = X0(4:6);

[a0,e0,~,~,~,~] = posvel2orbitalelements(r0,v0,p.mu);









%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Final BCs $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

rf = Xf(1:3);
vf = Xf(4:6);

[af,ef,~,~,~,~] = posvel2orbitalelements(rf,vf,p.mu);



%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Package Up $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


C_ineq = [e0 - p.e0;
          ef - p.ef;
          af - (p.af + 10000);
          (p.af - 10000) - af
          a0 - (p.a0 + 10000);
          (p.a0 - 10000) - a0];
          
C_eq = [];%[a0 - p.a0];
       % af - p.af];



end