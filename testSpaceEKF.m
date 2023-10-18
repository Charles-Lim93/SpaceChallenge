clear all; close all; clc;

dt = 0.05;
t = 0:dt:100;

Nsamples = length(t);

Xsaved = zeros(Nsamples, 6);
Zsaved = zeros(Nsamples, 3);

for k = 1:Nsamples
    [distance, yaw, pitch, target_x, target_y, target_z, target_vel_x, target_vel_y, target_vel_z] = GetData(dt);
    z = [distance, yaw, pitch]';
    disp("z @ test")
    disp(z)
    [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z] = SpaceEKF(z,dt);
    
    Xsaved(k, :) = [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z];
%     Xsaved(k, :) = [(relative_x+obs_x), (relative_y+obs_y), (relative_z+obs_z), vel_x, vel_y, vel_z];
%     Xsaved(k, :) = [(pos_x+obs_x), (pos_y+obs_y), (pos_z+obs_z), vel_x, vel_y, vel_z];
    Zsaved(k, :) = [distance, yaw, pitch]; 
    Tsaved(k, :) = [target_x, target_y, target_z, target_vel_x, target_vel_y, target_vel_z];
end

posXSaved = Xsaved(:,1);
posYSaved = Xsaved(:,2);
posZSaved = Xsaved(:,3);

velXSaved = Xsaved(:,4);
velYSaved = Xsaved(:,5);
velZSaved = Xsaved(:,6);

TrueXSaved = Tsaved(:,1);
TrueYSaved = Tsaved(:,2);
TrueZSaved = Tsaved(:,3);

TrueVelXSsaved = Tsaved(:,4);
TrueVelYSsaved = Tsaved(:,5);
TrueVelZSsaved = Tsaved(:,6);


t = 0:dt:Nsamples*dt - dt;

figure(1)
plot(t, velXSaved)
title("Velocity X")
hold on
plot(t, TrueVelXSsaved)
title("Velocity X")


figure(2)
plot(t, velYSaved)
hold on
plot(t, TrueVelYSsaved)
title("Velocity Y")

figure(3)
plot(t, velZSaved)
hold on
plot(t, TrueVelZSsaved)
title("Velocity Z")

figure(4)
plot(t, posXSaved, 'o-')
hold on
plot(t, TrueXSaved)
title("Position X")

figure(5)
plot(t, posYSaved, 'o-')
hold on
plot(t, TrueYSaved)
title("Position Y")

figure(6)
plot(t, posZSaved, 'o-')
hold on
plot(t, TrueZSaved)
title("Position Z")


