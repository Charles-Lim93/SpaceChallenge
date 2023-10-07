function [angle, pitch] = GetData(dt)

target_x = 100 + 0.5*dt + 1 * randn; % 가속도 텀?
target_y = 150 + 0.7*dt + 2 * randn; % 가속도 텀
target_z = 200 - 1.0*dt + 2 * randn; % + 가속도 텀 필요하지 않은지?

obs_x =150 - 1.0*dt; % + 가속도 텀
obs_y =120 - 1.0*dt; % + 가속도 텀
obs_z =100 - 1.0*dt; % + 가속도 텀

relative_x = (obs_x - target_x);
relative_y = (obs_y - target_y);
relative_z = obs_z - target_z;

disp(relative_x)
disp(relative_y)
disp(relative_z)

angle = atan2(relative_y, relative_x);
pitch = asind(-relative_z/sqrt(relative_x^2+relative_y^2+relative_z^2));

