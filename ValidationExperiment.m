% System:
l=1;
m=1;
g=9.81;
system=[l;m;g];
%Set the skill parameters for the seven skills
kps=[ 50; 60; 50; 120];
kds=[ 10; 8; 5; 4];

for skill = 1:length(kps)% Repteate the experiment of the diffrent skills
    %Set the initial values of the paramters for the Paramtrized loss

    % Set optimal skill parameters
    kp=kps(skill);
    kd=kds(skill);
    i=1000;
    OD1=zeros(i,2);% Set the Paramter Loss matrics
    OD2=zeros(i,2);% Set the General Loss matrics

    for n=1:i % Rpeate the experiment for 1000 trials

        [storing_list u_noise q_with_noise q_with_noise_dot] = Find_loss(system,kp,kd);

        % Locate the minimum loss optained by each risk metrics and return its value and location
        [loss3 ,location_of_Data4]= min(storing_list(:,3));
        [loss4 ,location_of_Data5]= min(storing_list(:,4));
        


        % append the parameter loss at the location of minimum risk loss to
        % the Parameter loss matrix
        OD1(n,1)=storing_list(location_of_Data4,4);
        OD1(n,2)=storing_list(location_of_Data5,4);
        

        % append the minimum risk loss to the General loss matrix
        OD2(n,1)=loss3;
        OD2(n,2)=loss4;


    end

    % Create a list to store the mean and the standard deviation of the
    % parameter loss assosiated with each risk metric
    list2={'Loss Metrics', 'mean Parameter Loss','Standard Deviation';
        'Obtimal Loss',sum(OD1(:,1))/i ,std(OD1(:,1));
        'Parameter Loss',sum(OD1(:,2))/i ,std(OD1(:,2))}

    % Create a list to store the mean of the risk loss assosiated with each risk metric
    list3={'Loss Metrics', 'Loss';
        'AEloss',sum(OD2(:,1))/i;
        'L_huber',sum(OD2(:,2))/i  }

    % Save the lists to an excel file
    filename = 'experminet3_PARAMETER_LOSS_new.xlsx';
    writecell(list2,filename,'Sheet',skill,'Range','D1')

    filename = 'experminet3_RISK_LOSS_new.xlsx';
    writecell(list3,filename,'Sheet',skill,'Range','D1')


end
function [storing_list u_noise q_with_noise q_with_noise_dot] = Find_loss(system,kp,kd)
% System:
l=system(1);
m=system(2);
g=system(3);
%Sampling
time=3;
ts=0.005;
t=time;
T=ts;

theta_opt=[kp;kd];
goal=pi/2;% set the goal
% Generate a noisy demonistration with noise with the mean of 1% of the states and actions value
[q_with_noise,q_with_noise_dot,u_noise] = get_demo(system, time, ts, theta_opt,0.01,goal);

%Get demonistration without the noise
[q,q_dot,u] = get_demo(system, time, ts, theta_opt,0,goal);


iteration=1000;
storing_list=zeros(iteration,4);% Create the sorting list as zero matrix

for i = 1:iteration  % Get 1000 line of data into sorting list

    %Data Selection, generate three data points in the range of number of samples
    while true
        N1=round((t/T)*rand(1));
        N2=round((t/T)*rand(1));
        if (N1~=0)&&(N2 ~= 0)
            break
        end
    end
    % Use the random data to select the data out of the noisy demo
    % create the feature vectors
    y_SD=[u_noise(N1);u_noise(N2)];
    phi_SD=[(goal-q_with_noise(N1)) (goal-q_with_noise(N2)) ;-q_with_noise_dot(N1) -q_with_noise_dot(N2)]';

    %Learning, Ridge regression
    theta=inv(phi_SD'*phi_SD+eye(2)*0.00000001)*phi_SD'*y_SD;

    %Data-Driven errors
    %Generate non-noisy demoinstration using the learned parameters
    [q_l,q_dot_l,u_l] = get_demo(system, time, ts, theta,0,goal);


    %Calculate the loss between the data from the learned
    %parameters and the nin-noisy demo using the skill parameters
    %Optimal data driven metric calculation
    
    A_l=((q_dot_l'*q_dot_l)*((goal*ones(length(q_l),1)-q_l)'*u_l))+(q_dot_l'*(goal*ones(length(q_l),1)-q_l)*(-q_dot_l'*u_l));
    D_l=(q_dot_l'*q_dot_l)*((goal*ones(length(q_l),1)-q_l)'*(goal*ones(length(q_l),1)-q_l))-(q_dot_l'*(goal*ones(length(q_l),1)-q_l))'*(q_dot_l'*(goal*ones(length(q_l),1)-q_l));
    B_l=((q_dot_l'*(goal*ones(length(q_l),1)-q_l))*((goal*ones(length(q_l),1)-q_l)'*u_l))+((goal*ones(length(q_l),1)-q_l)'*(goal*ones(length(q_l),1)-q_l))*(-q_dot_l'*u_l);
    A=((q_dot'*q_dot)*((goal*ones(length(q),1)-q)'*u))+(q_dot'*(goal*ones(length(q),1)-q)*(-q_dot'*u));
    D=(q_dot'*q_dot)*((goal*ones(length(q),1)-q)'*(goal*ones(length(q),1)-q))-(q_dot'*(goal*ones(length(q),1)-q))'*(q_dot'*(goal*ones(length(q),1)-q));
    B=((q_dot'*(goal*ones(length(q),1)-q))*((goal*ones(length(q),1)-q)'*u))+((goal*ones(length(q),1)-q)'*(goal*ones(length(q),1)-q))*(-q_dot'*u);
    P_data_l=[A_l/D_l; B_l/D_l];
    P_data=[A/D; B/D];
    data_loss = ((sum(abs(P_data-P_data_l))/(P_data(1)^2+P_data(2)^2)^0.5))*100;
    %Calculate the paramter loss
    loss=((sum(abs(theta-theta_opt))/(theta_opt(1)^2+theta_opt(2)^2)^0.5))*100;
    % Sotring Data
    storing_list(i,1)=N1;
    storing_list(i,2)=N2;
    storing_list(i,3)=data_loss;
    storing_list(i,4)=loss;



end


%Demo
    function [q_with_noise,q_with_noise_dot,u_noise] = get_demo(system, time, ts, theta,noize,goal)
        % System:
        l=system(1);
        m=system(2);
        g=system(3);
        num_steps = time / ts;

        %Initialisation
        qi=0;
        q_doti=0;
        q_with_noise=zeros(num_steps + 1, 1);
        q_with_noise_dot=zeros(num_steps + 1, 1);
        u_noise=zeros(num_steps + 1, 1);
        e=zeros(num_steps + 1, 1);
        q_with_noise(1)=qi;
        e(1)=goal-qi;
        q_with_noise_dot(1)=q_doti;
        u_noise(1)= theta'*[e(1) -q_with_noise_dot(1)]';

        for sys_iteration=1:num_steps

            e(sys_iteration)=goal-q_with_noise(sys_iteration);
            phi=[e(sys_iteration) -q_with_noise_dot(sys_iteration)]';
            u_noise(sys_iteration)=(theta'*phi)*(1+noize*randn(1));
            q_with_noise_dot(sys_iteration+1)=q_with_noise_dot(sys_iteration)+((u_noise(sys_iteration)/m*l.^2)-g*cos(q_with_noise(sys_iteration))/l)*ts;
            q_with_noise(sys_iteration+1)=q_with_noise(sys_iteration)+q_with_noise_dot(sys_iteration)*ts;
            q_with_noise_dot(sys_iteration)=q_with_noise_dot(sys_iteration)*(1-noize*rand(1));
            q_with_noise(sys_iteration)=q_with_noise(sys_iteration)*(1-noize*rand(1));
        end
    end
end