function [distance, angle, pitch, target_x, target_y, target_z, vel_x, vel_y, vel_z] = GetData(dt)

persistent initial_x initial_y initial_z initial_obs_x initial_obs_y initial_obs_z

if isempty(initial_x)
    initial_x = 100;
    initial_y = 100;
    initial_z = 100;
    initial_obs_x = 10;
    initial_obs_y = 10;
    initial_obs_z = 10;

end


vel_x = 5 + 1*dt;
vel_y = 10 + 1.5*dt;
vel_z = 20 + 1.5*dt;

vel_obs_x = -2;
vel_obs_y = 1;
vel_obs_z = 5;
target_x = initial_x + vel_x*dt ; % 가속도 텀?
target_y = initial_y + vel_y*dt ; % 가속도 텀
target_z = initial_z + vel_z*dt ; % + 가속도 텀 필요하지 않은지?

obs_x = initial_obs_x + vel_obs_x*dt; % 
obs_y = initial_obs_y + vel_obs_y*dt; % 
obs_z = initial_obs_z + vel_obs_z*dt; %

relative_x = target_x - obs_x;
relative_y = target_y - obs_y;
relative_z = target_z - obs_z;

v1 = 0.05*randn;
v2 = 0.05*randn;
v3 = 0.05*randn;

distance = sqrt(relative_x^2 + relative_y^2 + relative_z^2) + v1;

angle = (atan2(relative_y, relative_x)) + v2;

pitch = (asin(relative_z/distance)) + v3;


initial_x = target_x;
initial_y = target_y;
initial_z = target_z;