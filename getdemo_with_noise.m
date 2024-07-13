
function [q1 q2 q1dot q2dot int_e1 int_e2 u1 u2 ] = getdemo_with_noise(time,ts,phi1,phi2,system,q1f,q2f,noise)

kp1 = phi1(1);
kp2 = phi2(1);
ki1 = phi1(2);
ki2 = phi2(2);
kd1 = phi1(3);
kd2 = phi2(3);

m1 = system(1);l1 = system(2);m2 = system(3);l2 = system(4);g = system(5);

% Preallocate arrays with one extra element
num_steps = time / ts;
q1 = zeros(num_steps + 1, 1);
q2 = zeros(num_steps + 1, 1);
q1dot = zeros(num_steps + 1, 1);
q2dot = zeros(num_steps + 1, 1);
q1ddot = zeros(num_steps + 1, 1);
q2ddot = zeros(num_steps + 1, 1);
e1 = zeros(num_steps + 1, 1);
e2 = zeros(num_steps + 1, 1);
int_e1= zeros(num_steps + 1, 1);
int_e2= zeros(num_steps + 1, 1);
u1= zeros(num_steps + 1, 1);
u2= zeros(num_steps + 1, 1);
% Initial conditions
q1(1) = 0;
q2(1) = 0;
% Final desired positions
q1f = q1f;
q2f = q2f;

e1(1)=q1f - q1(1);
e2(1)=q2f - q2(1);


for i = 1:num_steps
    % Compute the control torques
    e1(i)=q1f - q1(i);
    e2(i)=q2f - q2(i);
    int_e1(i)=sum(e1(1:i))*ts;
    int_e2(i)=sum(e2(1:i))*ts;
    T1 = kp1 * (e1(i)) - kd1 * q1dot(i) + ki1 * int_e1(i);
    T2 = kp2 * (e2(i)) - kd2 * q2dot(i) + ki2 * int_e2(i);
    u1(i)=T1*(1+noise*rand(1));
    u2(i)=T2*(1+noise*rand(1));
    T = [u1(i); u2(i)];

    % Dynamics matrices
    M = [m1*(l1^2) + m2*((l1^2) + 2*l1*l2*cos(q2(i)) + l2^2), m2*(l1*l2*cos(q2(i)) + l2^2);
        m2*(l1*l2*cos(q2(i)) + l2^2), m2*l2^2];
    C = [-m2*l1*l2*sin(q2(i))*(2*q1dot(i)*q2dot(i) + q2dot(i)^2);
        m2*l1*l2*q1dot(i)^2*sin(q2(i))];
    G = [(m1+m2)*l1*g*cos(q1(i)) + m2*g*l2*cos(q1(i) + q2(i));
        m2*g*l2*cos(q1(i) + q2(i))];
    if (det(M)~=0 && det(M)> -1000000000 && det(M) < 10000000000)

    % Acceleration
    a_r = M \ (T - C - G);
    q1ddot(i + 1) = a_r(1);
    q2ddot(i + 1) = a_r(2);

    % Velocity
    q1dot(i + 1) = q1dot(i) + q1ddot(i + 1) * ts;
    q2dot(i + 1) = q2dot(i) + q2ddot(i + 1) * ts;

    % Position
    q1(i + 1) = q1(i) + q1dot(i + 1) * ts;
    q2(i + 1) = q2(i) + q2dot(i + 1) * ts;
    q1(i + 1) = q1(i + 1) * (1+noise*rand(1));
    q2(i + 1) = q2(i + 1) * (1+noise*rand(1));
    end 

end
end