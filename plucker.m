
lambda = 1;

p11 = [lambda 12 10 10]';
p12 = [lambda 12 12 12]';
p1 = [284 240 1];

p21 = [lambda 20 20 20]';
p22 = [lambda 25 25 25]';
p2 = [350 230 1];

p31 = [lambda 30 30 30]';
p32 = [lambda 35 35 35]';
p3 = [360 220 1];

p41 = [lambda 40 40 40]';
p42 = [lambda 45 45 45]';
p4 = [370 265 1];

p51 = [lambda 50 50 50]';
p52 = [lambda 55 55 55]';
p5 = [370 280 1];

p61 = [lambda 60 60 60]';
p62 = [lambda 65 65 65]';
p6 = [380 220 1];

M = [1 0 0 0;
     0 1 0 0;
     0 0 1/2 0]; % 3 X 4, Camera Projection Matrix

L = plucker_calculation(p11, p12); % Plucker Matrix

l1 = cross(M*p11, M*p12);
l2 = cross(M*p21, M*p22);
l3 = cross(M*p31, M*p32);
l4 = cross(M*p41, M*p42);
l5 = cross(M*p51, M*p52);
l6 = cross(M*p61, M*p62);


M_tilt = [
                             plucker_calculation(M(2,:), M(3,:));
                             plucker_calculation(M(3,:), M(1,:));
                             plucker_calculation(M(1,:), M(2,:))
                             ];

a = p1 * M_tilt;
b = p2 * M_tilt;
c = p3 * M_tilt;
d = p4 * M_tilt;
e = p5 * M_tilt;
f = p6 * M_tilt;

K = [a; b; c; d; e; f];
B = [0; 0; 0; 0; 0; 0];
disp(inv(K))
disp("ANS")
X = inv(K) * B;
% syms dx dy dz mx my mz
% 
% 
% dx, dy, dz, mx, my, mz = linsolve(K,B)
% 
% 
 