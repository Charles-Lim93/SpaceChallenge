

dt = 0.05;
t = 0:dt:20;

Nsamples = length(t);

Xsaved = zeros(Nsamples, 6);
Zsaved = zeros(Nsamples, 2);

for k = 1:Nsamples
    [yaw, pitch] = GetData(dt);
    z = [yaw, pitch];

    [vel_x, vel_y, vel_z, pos_x, pos_y, pos_z] = SpaceEKF(z,dt);
    
    Xsaved(k, :) = [vel_x vel_y vel_z pos_x pos_y pos_z];
    Zsaved(k, :) = [yaw, pitch]; 
end

velXSaved = Xsaved(:,1);
velYSaved = Xsaved(:,2);
velZSaved = Xsaved(:,3);
posXSaved = Xsaved(:,4);
posYSaved = Xsaved(:,5);
posZSaved = Xsaved(:,6);

t = 0:dt:Nsamples*dt - dt;

figure(1)
plot(t, velXSaved)
title("Velocity X")
figure(2)
plot(t, velYSaved)
title("Velocity Y")

figure(3)
plot(t, velZSaved)
title("Velocity Z")

figure(4)
plot(t, posXSaved)
title("Position X")

figure(5)
plot(t, posYSaved)
title("Position Y")

figure(6)
plot(t, posZSaved)
title("Position Z")


