function [v0_v,v_v] = Universal_Lambert(r0_v,r_v,TOF,tm,mu)

    r0=norm(r0_v);
    r=norm(r_v);
    
    cosdnu=dot(r0_v,r_v)/(r*r0);

    A=tm*sqrt(r*r0*(1+cosdnu));

    if A == 0
        error('Universal_Lambert: A is zero')
    end

    psi_n=0;
    c2=0.5;
    c3=1/6;

    psi_up=4*pi^2;
    psi_low = -4*pi^2;
    dtn = TOF*100;
    count = 0;

    while abs(TOF-dtn) > 1e-6 && count < 1000
        yn=r0 + r + (A*(psi_n*c3 - 1))/sqrt(c2);

        % if (A > 0) && (yn < 0)
        %     %need to find lower bound such that y is positive
        %     psi_low = fzero(@(x) posypsibnd(x,A,r0,r),0);
        %     psin = (psi_low + psi_up)/2;
        % 
        % 
        %     if psin >=0
        %         c2 = (1 -cos(sqrt(psin)))/psin;
        %         c3 = (sqrt(psin) - sin(sqrt(psin)))/(sqrt(psin)^3);
        %     else
        %         c2 = (1 -cosh(sqrt(-psin)))/psin;
        %         c3 = (-sqrt(-psin) + sin(sqrt(-psin)))/(sqrt(-psin)^3);
        %     end
        % 
        %     yn = r0 + r + A*(psin*c3 - 1)/sqrt(c2);
        % end

        if yn < 0 
            psi_low = psi_low+1;
        end

        chi_n = sqrt(yn/c2);
        dtn = (c3*chi_n^3 + A*sqrt(yn))/sqrt(mu);

        if dtn<TOF
            psi_low = psi_n;
        else
            psi_up = psi_n;
        end

        psi_np1 = 0.5*(psi_up + psi_low);

        if psi_np1 >=0
            c2 = (1 -cos(sqrt(psi_np1)))/psi_np1;
            c3 = (sqrt(psi_np1) - sin(sqrt(psi_np1)))/(sqrt(psi_np1)^3);
        else
            c2 = (1 -cosh(sqrt(-psi_np1)))/psi_np1;
            c3 = (-sqrt(-psi_np1) + sin(sqrt(-psi_np1)))/(sqrt(-psi_np1)^3);
        end
        psi_n = psi_np1;
        count = count + 1;



    end

    f=1-(yn/r0);
    g_dot = 1-(yn/r);
    g=A*sqrt(yn/mu);

    v0_v = (r_v - f*r0_v)/g;
    v_v = (g_dot*r_v - r0_v)/g;


%helper function to calculate c2/c3
     function [c2,c3] = c2c3(psi)
         
         c2 = zeros(size(psi));
         c3 = zeros(size(psi));
         
         zeropsi = psi == 0;
         pospsi = psi > 0;
         negpsi = psi < 0;
         
         c2(zeropsi) = 1/2;
         c3(zeropsi) = 1/6;
         
         c2(pospsi) = (1 - cos(sqrt(psi(pospsi))))./psi(pospsi);
         c3(pospsi) = (sqrt(psi(pospsi)) - sin(sqrt(psi(pospsi))))./...
             sqrt(psi(pospsi).^3);
         
         c2(negpsi) = (1 - cosh(sqrt(-psi(negpsi))))./psi(negpsi);
         c3(negpsi) = (sinh(sqrt(-psi(negpsi))) - sqrt(-psi(negpsi)))./...
             sqrt(-psi(negpsi).^3);
         
     end
 
     %helper function for positive y bound finding
     function del = posypsibnd(psi,A,r0,r)
         
         [c2,c3] = c2c3(psi);
         del = (psi*c3 - 1)/sqrt(c2) + (r0+r)/A;
         
     end


end
