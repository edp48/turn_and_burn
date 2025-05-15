function vec_out = rot_z(vec_in,theta_in)
%This function takes in a vector and rotates it by a value of positive 
%theta (or negative if the sign of theta is negative) as given by the
%standard DCM for a z_axis rotation
%
%vec_in must be a 3x1 col vector
%theta is in degrees

sz = size(vec_in);
if sz ~= [3,1]
    error("Wrong Size Vector")
end 

theta = theta_in; %deg2rad(theta_in);

R_z = [cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];

vec_out = R_z*vec_in;


end 