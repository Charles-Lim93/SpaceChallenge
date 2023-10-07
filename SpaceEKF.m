function [vel_x, vel_y, vel_z, pos_x, pos_y, pos_z] = SpaceEKF(z, dt)


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

    R = 100 * eye(6);

    x = [2 2 3 50 100 100]'; % Initial Value
    P = 10 * eye(2);

    firstRun = 1;
end

H = Hjacob(x);

xp = A*x; % 6 X 1
Pp = A*P*A' + Q; % 6 X 6
K = Pp*H'*inv(H*Pp*H' + R); 

x = xp + K*(z - hx(xp));
P = Pp - K*H*Pp;

vel_x = x(1);
vel_y = x(2);
vel_z = x(3);
pos_x = x(4);
pos_y = x(5);
pos_z = x(6);

%------------------------
% Observation Value

function zp = hx(xhat)

x1_hat = xhat(1);
x2_hat = xhat(2);
x3_hat = xhat(3);
zp = [atan2(x2_hat,x1_hat) asind(x3_hat/sqrt(x1_hat^2 + x2_hat^2 + x3_hat^2))]';
%-------------
% JACOBIAN
function H = Hjacob(xp)
    x1 = xp(4); % position x
    x2 = xp(5); % position y
    x3 = xp(6); % position z

    H = zeros(2, 6);
    distance = sqrt(x1^2 + x2^2 + x3^2);

    H(1,4:6) = [-x2/(x1^2 + x2^2), x1/(x1^2 + x2^2), 0];
    H(2,4:6) = [-(x1*x3)/(sqrt(x1^2+x2^2)*distance), -(x2*x3)/(sqrt(x1^2+x2^2)*distance), ...
    sqrt(x1^2+x2^2)/distance];


