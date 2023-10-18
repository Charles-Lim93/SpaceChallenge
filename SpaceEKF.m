function [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z] = SpaceEKF(z, dt)


persistent A Q R
persistent x P
persistent firstRun

if isempty(firstRun)
    A = eye(6) + dt * [0 0 0 1 0 0;
                        0 0 0 0 1 0;
                        0 0 0 0 0 1;
                        0 0 0 0 0 0;
                        0 0 0 0 0 0;
                        0 0 0 0 0 0];

    B = zeros(6,3) + dt * [0 0 0;
                            0 0 0;
                            0 0 0;
                            -1 0 0;
                            0 -1 0;
                            0 0 -1];

    Q = [0.001 0 0 0 0 0;
        0 0.001 0 0 0 0;
        0 0 0.001 0 0 0;
        0 0 0 0.001 0 0;
        0 0 0 0 0.001 0;
        0 0 0 0 0 0.001];

    R = [2 0 0;
         0 0.5 0;
         0 0 0.5];

    x = [100 100 100 5 10 15]'; % Initial Value
    P = 15 * eye(6);

    firstRun = 1;
end
disp("X in SPACE EKF")
disp(x)

H = Hjacob(x);
disp("H")
disp(H)
%xp = A*x + B*u; % 6 X 1  Accelerating
xp = A*x; % 6 X 1

Q = A*Q*A*dt;
Pp = A*P*A' + Q; % 6 X 6
K = Pp * H'/ (H * Pp * H' + R); 
% disp("A")
% disp(A)
% disp("x")
% disp(x)
% disp("xp")
% disp(xp)
% disp("Pp")
% disp(Pp)
% disp("K")
% disp(K)
% disp("z")
% disp(z)
% disp("hx(xp)")
x = xp + K*(z - hx(xp));
disp(hx(xp))

disp("new x")
disp(x)
% P = Pp - K*H*Pp;

% P = Pp - K*H*Pp;
P = (eye(6)-K*H)*Pp*(eye(6)-K*H)' + K*R*K';
disp("new P")
disp(P)
disp(size(x));  % Should display  6     1
disp(size(xp));  % Should display  6     1
disp(size(K));   % Should display  6     3
disp(size(z));        % Should display  3     1
disp(size(hx(xp)));   % Should display  3     1
disp(size(H))  % Should display  3     6
pos_x = x(1);
pos_y = x(2);
pos_z = x(3);
vel_x = x(4);
vel_y = x(5);
vel_z = x(6);

%------------------------
% Observation Value

function zp = hx(xhat)
disp("xhat in hx")
disp(xhat)
x1_hat = xhat(1); % X
x2_hat = xhat(2); % Y
x3_hat = xhat(3); % Z
distance = sqrt(x1_hat^2 + x2_hat^2 + x3_hat^2);
angle = atan2(x2_hat, x1_hat);
pitch = asin(x3_hat / distance);
% angle = rad2deg(angle);
zp = [distance, (angle),  (pitch)]';

%-------------
% JACOBIAN
function H = Hjacob(xp)
    disp("xp in Jacobian")
    disp(xp)
    x1 = (xp(1)); % position x
    x2 = (xp(2)); % position y
    x3 = (xp(3)); % position z

    H = zeros(3, 6);
    d = sqrt(x1^2 + x2^2 + x3^2);

    H(1, 1:3) = [x1/d, x2/d, x3/d];
    H(2, 1:3) = [-x2/(sqrt(x1^2+x2^2)), x1/(sqrt(x1^2+x2^2)), 0];
    H(3, 1:3) = [(-x1*x3)/(sqrt(x1^2+x2^2+x3^2)*d^2) (-x2*x3)/(sqrt(x1^2+x2^2+x3^2)*d^2) (x1^2+x2^2)/(sqrt(x1^2+x2^2+x3^2)*d^2)];


