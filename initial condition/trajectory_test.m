clc;clear;
close all;

%% Earth

a = 6378.137;
flattening = 1/298.257223563;
b = a*(1-flattening);

[X,Y,Z] = ellipsoid(0,0,0,a,a,b);

figure
earth = surf(X,Y,Z,'FaceAlpha',0);
earth.EdgeColor = [0.7,0.7,0.7];
earth.LineStyle = ":";
hold on
axis equal

k = 0 : 0.1*pi : 2*pi;
equaterX = 8000*cos(k);
equaterY = 8000*sin(k);

equator = patch(equaterX,equaterY,'c','FaceAlpha',0.2);
equator.EdgeColor = [1 1 1];

quiver3(0,0,0,1,0,0,7500);
quiver3(0,0,0,0,1,0,7500);
quiver3(0,0,0,0,0,1,7500);

text(8000, 0, 0, 'X_{ECI}')
text(0, 8000, 0, 'Y_{ECI}')
text(0, 0, 8000, 'Z_{ECI}')

mu_e = 3.986004418*10^14;  % [m^3/s^2]
mu_e = mu_e/1e9;


%% Earth rate // rad/sec
earth_rate = 7.29211585537707E-05;

%% Missile apogee, flight time

missile_h = 97;
flight_time = 4.8147;

%% Launch Point 

r_launch_ecef = lla2ecef([39.7437 127.4732 0], 'WGS84');
r_launch_ecef = r_launch_ecef/1e3;

r_launch_eci = r_launch_ecef';
plot3(r_launch_ecef(1),r_launch_ecef(2),r_launch_ecef(3),'o','LineWidth',3);
% text(r_launch(1)+500,r_launch(2)+500,r_launch(3)+500,'Launch')

%% Land Point

r_land_ecef = lla2ecef([36.37247 127.3578 0], 'WGS84');
r_land_ecef = r_land_ecef/1e3;

r_land_eci = ecef2eci(flight_time*60, earth_rate)*r_land_ecef';
% plot3(r_land_eci(1),r_land_eci(2),r_land_eci(3),'*');
% text(r_land(1)+500,r_land(2)+500,r_land(3)+500,'Land')

%% Time // 5min

dt = 0.1;
t = 0 : dt : flight_time*60;
% t= 0 : dt : 90;

%% Missile Trajectory Reconstruction 

% apogee point
mid_eci = (r_launch_eci + r_land_eci)/2;
mid_ecef = ecef2eci(flight_time*60/2, earth_rate)'*mid_eci;
mid_lla = ecef2lla(1e3*mid_ecef');
mid_lla(3) = missile_h*1e3; 
ap_ecef = lla2ecef(mid_lla)';
ap_ecef = ap_ecef/1e3;
ap_eci = ecef2eci(flight_time*60/2, earth_rate)*ap_ecef;

plot3(ap_eci(1),ap_eci(2),ap_eci(3),'*','LineWidth',3);

A_z_vec = cross(r_launch_eci,ap_eci)/norm(cross(r_launch_eci,ap_eci));
A_x_vec = -ap_eci/norm(ap_eci);
A_y_vec = cross(A_z_vec, A_x_vec);

% quiver3(0,0,0,A_z_vec(1),A_z_vec(2),A_z_vec(3),7500);
% quiver3(0,0,0,A_x_vec(1),A_x_vec(2),A_x_vec(3),7500);
% quiver3(0,0,0,A_y_vec(1),A_y_vec(2),A_y_vec(3),7500);

dcm_eci2plane = [A_x_vec A_y_vec A_z_vec]';

launch_plane = dcm_eci2plane*r_launch_eci;
ap_plane = dcm_eci2plane*ap_eci;

ra = abs(ap_plane(1));
nu = pi - acos(dot(ap_plane,launch_plane)/(norm(ap_plane)*norm(launch_plane)));
e_missile = (ra-norm(launch_plane))/(ra+norm(launch_plane)*cos(nu));
a_missile = ra/(1+e_missile);

v_ap_plane = [0; -sqrt(2*mu_e*(1/ra - 1/(2*a_missile))); 0];
v_launch_plane_norm = sqrt(2*mu_e*(1/norm(launch_plane) - 1/(2*a_missile)));

slope = -(1-e_missile^2)*(launch_plane(1)+a_missile*e_missile)/launch_plane(2);
v_launch_plane = [-v_launch_plane_norm/sqrt(1+slope^2) ;-(slope)*v_launch_plane_norm/sqrt(1+slope^2);0];

v_launch_eci = dcm_eci2plane'*v_launch_plane;

%% RK4
x_missile(:,1) = [r_launch_eci ; v_launch_eci];
x_missile_ecef(:,1) = x_missile(1:3,1);
u = 0;
for i = 1 : size(t,2)-1
    x_missile(:,i+1) = RK4(x_missile(:,i),dt,u,mu_e);
    time = i*dt;
    
    x_missile_ecef(:,i+1) = ecef2eci(time,earth_rate)'*x_missile(1:3,i+1);
    x_missile_lla = ecef2lla(1e3*x_missile_ecef(:,i+1)');

    % if x_missile_lla(3) <= 0
    %     % disp(i/(10*60));
    %     break
    % end
end

plot3(x_missile(1,:),x_missile(2,:),x_missile(3,:),'LineWidth',2);
plot3(x_missile(1,size(x_missile,2)),x_missile(2,size(x_missile,2)),x_missile(3,size(x_missile,2)),'x','LineWidth',3)


r_launch_lla = ecef2lla(1e3*r_launch_eci');
r_ap_lla = ecef2lla(1e3*ap_eci');
x_missile_lla = ecef2lla(1e3*x_missile_ecef');

%% Satellite // VICTUS NOX
tle = loadTle(pwd, 'Victus_Nox.tle');

for i = 0 : 5
    [R(i+1,:), V(i+1,:)] = sgp4(tle, 58002+i);
    R_lla = ecef2lla(1e3*R(i+1,:));
    % disp(R_lla(3))

end

plot3(R(:,1),R(:,2),R(:,3),'*')
% [R1, V1] = sgp4(tle, 25);
% plot3(R1(1),R1(2),R1(3),'*')

%% RK4 

x_sat = [R(1,:)';V(1,:)'];
x_sat_ecef(:,1) = x_sat(1:3,1);
x_sat_lla(1,:) = ecef2lla(1e3*x_sat_ecef(:,1)');
u = 0;

for i = 1 : size(t,2)-1
    x_sat(:,i+1) = RK4(x_sat(:,i),dt,u,mu_e);
    time = i*dt;
    
    x_sat_ecef(:,i+1) = ecef2eci(time,earth_rate)'*x_sat(1:3,i+1);
    x_sat_lla(i+1,:) = ecef2lla(1e3*x_sat_ecef(:,i+1)');
end

plot3(x_sat(1,:),x_sat(2,:),x_sat(3,:),'LineWidth',2);

%% camera center // seoul 37.566535 126.9779692
% r_launch_ecef = lla2ecef([39.7437 127.4732 0], 'WGS84');
r_cam_center_ecef = lla2ecef([39.7437 127.4732 0], 'WGS84');

for i = 0 : size(t,2)-1
    r_cam_center_eci(:,i+1) = 1e-3*ecef2eci(i*dt, earth_rate)*r_cam_center_ecef';
end

for i = 0 : size(t,2)-1
    cam_missile_eci(:,i+1) = x_missile(1:3,i+1) - x_sat(1:3,i+1);
    cam_center_eci_norm(:,i+1) = (r_cam_center_eci(:,i+1) - x_sat(1:3,i+1))/norm(r_cam_center_eci(:,i+1) - x_sat(1:3,i+1)); 
    cam_missile_eci_norm(:,i+1) = (x_missile(1:3,i+1) - x_sat(1:3,i+1))/norm(x_missile(1:3,i+1) - x_sat(1:3,i+1));
    % angle(i+1)=acosd(dot(cam_missile_eci(:,i+1),cam_center_eci(:,i+1)));
end

% figure;
% plot(t,angle);

%eci 2 vvlh 
for i = 0 : size(t,2)-1
    dcm_eci2vvlh(:,:,i+1) = eci2vvlh(x_sat(1:3,i+1),x_sat(4:6,i+1));
    cam_center_vvlh(:,i+1) = dcm_eci2vvlh(:,:,i+1)*cam_center_eci_norm(:,i+1);
    dcm_vvlh2cam(:,:,i+1) = vvlh2cam(cam_center_vvlh(:,i+1));
    dcm_eci2cam(:,:,i+1) = dcm_vvlh2cam(:,:,i+1)*dcm_eci2vvlh(:,:,i+1);
    cam_missile_cam(:,i+1) = dcm_eci2cam(:,:,i+1)*cam_missile_eci(:,i+1);
end

for i = 0 : size(t,2)-1
    dcm_eci2cam(:,:,i+1)*cam_missile_eci_norm(:,i+1);
end

%camera matrix
f = 100; %[mm]
f = f*1e-3; %[m]
Cx = 0.5*5328;
Cy = 0.5*4608;
m_pixel = 1/(2.74*1e-6); %[um]

M_cam = [m_pixel*f  0 Cx 0; 
         0  m_pixel*f Cy 0;
         0    0       1  0];

for i = 0 : size(t,2)-1
    missile_img(:,i+1) = M_cam*[1e3*cam_missile_cam(:,i+1);1];
    missile_img(:,i+1) =  round(missile_img(:,i+1)/missile_img(3,i+1));
end

figure
plot(missile_img(1,:),missile_img(2,:));
hold on
plot(Cx,Cy,'x');
xlim([0 5328]);
ylim([0 4608]);

%%save file

save('Ground_Truth1.mat','t','r_launch_eci','x_missile','x_sat','dcm_eci2cam','flight_time','missile_img','-v7.3');

%% Geoplot

uif = uifigure;
g = geoglobe(uif);

% geoplot3(g,r_launch_lla(1),r_launch_lla(2),r_launch_lla(3),'o','LineWidth',5,'Color','y');
% hold(g,'on')
% geoplot3(g,r_ap_lla(1),r_ap_lla(2),r_ap_lla(3),'o','LineWidth',5,'Color','y');
geoplot3(g,x_missile_lla(:,1)',x_missile_lla(:,2)',x_missile_lla(:,3)','LineWidth',2,'Color','y');
% geoplot3(g,x_lla(1,size(x_lla,2)-1),x_lla(2,size(x_lla,2)-1),x_lla(3,size(x_lla,2)-1),'o','LineWidth',5,'Color','y')

hold(g,'on')
geoplot3(g,x_sat_lla(:,1)',x_sat_lla(:,2)',x_sat_lla(:,3)','LineWidth',2,'Color','c');
%% Function

%% ecef2eci 
function dcm = ecef2eci(t,rate) % t in sec
    theta = rate * t;
    
    c = cos(theta);
    s = sin(theta);

    dcm = [c -s 0 ;
           s c 0 ;
           0 0 1];
end

%% RK4 
% input : r = state / h = time step / u = input_accel = 0 
%         / mu = gravitational constant
% output : r_out = result of RK4
% RK4 function
function r_out = RK4(r,h,u,mu)
    
    k1 = dyn_eqn(r,u,mu);
    k2 = dyn_eqn(r+0.5*k1*h,u,mu);
    k3 = dyn_eqn(r+0.5*k2*h,u,mu);
    k4 = dyn_eqn(r+k3*h,u,mu);
    
    r_out = r + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end

%% Dynamic Eq
% input : r = state / u = input_accel = 0 / mu = gravitational constant
% output : r_dot 
% differential eq represting dynamic of system
function r_dot = dyn_eqn(r,u, mu)

    R = [r(1) ; r(2) ; r(3)];
    R_dot = [r(4) ; r(5) ; r(6)];

    R_dotdot = -(mu/norm(R)^3)*R + u;

    r_dot = [R_dot ; R_dotdot];
end

%% ECI to vvlh

function dcm = eci2vvlh(r,v)
    z = -r/norm(r);
    y = cross(z,v)/norm(cross(z,v));
    x = cross(y,z);

    dcm = [x y z]';
end

%% vvlh to cam // z = cam center / y = inv velocity direction

function dcm = vvlh2cam(cam)

    z = cam/norm(cam);
    x = -cross(z,[-1;0;0])/norm(cross(z,[-1;0;0]));
    y = cross(z,x);

    dcm = [x y z]';
end