% Define your A and C matrices (Jacobian of system dynamics and observation model)

A = [1 0 0 1 0 0;
    0 1 0 0 1 0;
    0 0 1 0 0 1;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];
C = [sqrt(3)/3 sqrt(3)/3 sqrt(3)/3 0 0 0;
    -0.005 0.005 0 0 0 0;
    -sqrt(3)/9 -sqrt(3)/9 2*sqrt(3)/9 0 0 0];...; % Replace with your C matrix

% Number of states
n = size(A, 1);

% Calculate the observability matrix
O = C;
for i = 1:n-1
    O = [O; C * A^i];
end
disp(rank(O))
% Check observability
if rank(O) == n
    disp('The system is observable.');
else
    disp('The system is not observable.');
end

ob = obsv(A,C)

ubobsv = length(A) - rank(ob)