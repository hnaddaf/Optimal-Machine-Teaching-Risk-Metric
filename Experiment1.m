% System:
l=1;
m=1;
g=9.81;
system=[l;m;g];
%Set the skill parameters for the seven skills
kps=[ 178; 400; 94; 400; 170; 81; 40];
kds=[ 24; 50; 22; 50; 20; 12; 12];
kis=[ 20; 4000; 20; 800; 220; 320; 120];



for skill = 1:length(kps)% Repteate the experiment of the diffrent skills
    %Set the initial values of the paramters for the Paramtrized loss

    % Set optimal skill parameters
    kp=kps(skill);
    kd=kds(skill);
    ki=kis(skill);
    i=1000;
    OD1=zeros(i,11);% Set the Paramter Loss matrics
    OD2=zeros(i,11);% Set the General Loss matrics

    for n=1:i % Rpeate the experiment for 1000 trials

        [storing_list u_noise q_with_noise integral_e q_with_noise_dot] = Find_loss(system,kp,kd,ki);

        % Locate the minimum loss optained by each risk metrics and return its value and location
        [loss4 ,location_of_Data4]= min(storing_list(:,4));
        [loss5 ,location_of_Data5]= min(storing_list(:,5));
        [loss6 ,location_of_Data6]= min(storing_list(:,6));
        [loss7 ,location_of_Data7]= min(storing_list(:,7));
        [loss8 ,location_of_Data8]= min(storing_list(:,8));
        [loss9 ,location_of_Data9]= min(storing_list(:,9));
        [loss10 ,location_of_Data10]= min(storing_list(:,10));
        [loss11 ,location_of_Data11]= min(storing_list(:,11));
        [loss12 ,location_of_Data12]= min(storing_list(:,12));
        [loss13 ,location_of_Data13]= min(storing_list(:,13));


        % append the parameter loss at the location of minimum risk loss to
        % the Parameter loss matrix
        OD1(n,1)=storing_list(location_of_Data4,14);
        OD1(n,2)=storing_list(location_of_Data5,14);
        OD1(n,3)=storing_list(location_of_Data6,14);
        OD1(n,4)=storing_list(location_of_Data7,14);
        OD1(n,5)=storing_list(location_of_Data8,14);
        OD1(n,6)=storing_list(location_of_Data9,14);
        OD1(n,7)=storing_list(location_of_Data10,14);
        OD1(n,8)=storing_list(location_of_Data11,14);
        OD1(n,9)=storing_list(location_of_Data12,14);
        OD1(n,10)=storing_list(location_of_Data13,14);
        OD1(n,11)=min(storing_list(:,14));

        % append the minimum risk loss to the General loss matrix
        OD2(n,1)=loss4;
        OD2(n,2)=loss5;
        OD2(n,3)=loss6;
        OD2(n,4)=loss7;
        OD2(n,5)=loss8;
        OD2(n,6)=loss9;
        OD2(n,7)=loss10;
        OD2(n,8)=loss11;
        OD2(n,9)=loss12;
        OD2(n,10)=loss13;
        OD2(n,11)=min(storing_list(:,14));


    end

    % Create a list to store the mean and the standard deviation of the
    % parameter loss assosiated with each risk metric
    list2={'Loss Metrics', 'mean Parameter Loss','Standard Deviation';
        'AEloss',sum(OD1(:,1))/i ,std(OD1(:,1));
        'L_huber',sum(OD1(:,2))/i ,std(OD1(:,2)) ;
        'MSLE',sum(OD1(:,3))/i ,std(OD1(:,3));
        'NMSE',sum(OD1(:,4))/i ,std(OD1(:,4));
        'RMSE',sum(OD1(:,5))/i ,std(OD1(:,5));
        'RMSE_u',sum(OD1(:,6))/i ,std(OD1(:,6)) ;
        'NMSE_u',sum(OD1(:,7))/i ,std(OD1(:,7)) ;
        'AError',sum(OD1(:,8))/i ,std(OD1(:,8));
        'MSLE_u',sum(OD1(:,9))/i ,std(OD1(:,9));
        'Leaned_loss',sum(OD1(:,10))/i ,std(OD1(:,10));
        'Parameter_loss',sum(OD1(:,11))/i ,std(OD1(:,11))}

    % Create a list to store the mean of the risk loss assosiated with each risk metric
    list3={'Loss Metrics', 'Loss';
        'AEloss',sum(OD2(:,1))/i;
        'L_huber',sum(OD2(:,2))/i  ;
        'MSLE',sum(OD2(:,3))/i ;
        'NMSE',sum(OD2(:,4))/i ;
        'RMSE',sum(OD2(:,5))/i ;
        'RMSE_u',sum(OD2(:,6))/i ;
        'NMSE_u',sum(OD2(:,7))/i ;
        'AError',sum(OD2(:,8))/i;
        'MSLE_u',sum(OD2(:,9))/i;
        'Leaned_loss',sum(OD2(:,10))/i ;
        'Parameter_loss',sum(OD2(:,11))/i}

    % Save the lists to an excel file
    filename = 'experminet1_PARAMETER_LOSS_new.xlsx';
    writecell(list2,filename,'Sheet',skill,'Range','D1')

    filename = 'experminet1_RISK_LOSS_new.xlsx';
    writecell(list3,filename,'Sheet',skill,'Range','D1')

    data_bar = [sum(OD1(:,1))/i, sum(OD1(:,2))/i, sum(OD1(:,3))/i, sum(OD1(:,4))/i, sum(OD1(:,5))/i, sum(OD1(:,6))/i, sum(OD1(:,7))/i, sum(OD1(:,8))/i, sum(OD1(:,9))/i, sum(OD1(:,10))/i, sum(OD1(:,11))/i];


    labels = {'AEloss','L_huber','MSLE','NMSE','RMSE','RMSE_u','NMSE_u','AError','MSLE_u','Leaned_loss', 'paramter loss'};


end
function [storing_list u_noise q_with_noise integral_e q_with_noise_dot] = Find_loss(system,kp,kd,ki)
% System:
l=system(1);
m=system(2);
g=system(3);
%Sampling
time=3;
ts=0.005;
t=time;
T=ts;

theta_opt=[kp;ki;kd];
goal=pi/2;% set the goal
% Generate a noisy demonistration with noise with the mean of 1% of the states and actions value
[q_with_noise,q_with_noise_dot,u_noise,integral_e] = get_demo(system, time, ts, theta_opt,0.01,goal);

%Get demonistration without the noise
[q,q_dot,u] = get_demo(system, time, ts, theta_opt,0,goal);


iteration=1000;
storing_list=zeros(iteration,14);% Create the sorting list as zero matrix

for i = 1:iteration  % Get 1000 line of data into sorting list

    %Data Selection, generate three data points in the range of number of samples
    while true
        N1=round((t/T)*rand(1));
        N2=round((t/T)*rand(1));
        N3=round((t/T)*rand(1));
        if (N1~=0)&&(N2 ~= 0)&&(N3~=0)
            break
        end
    end
    % Use the random data to select the data out of the noisy demo
    % create the feature vectors
    y_SD=[u_noise(N1);u_noise(N2);u_noise(N3)];
    phi_SD=[(goal-q_with_noise(N1)) (goal-q_with_noise(N2)) (goal-q_with_noise(N3));integral_e(N1) integral_e(N2) integral_e(N3) ;-q_with_noise_dot(N1) -q_with_noise_dot(N2)  -q_with_noise_dot(N3)]';

    %Learning, Ridge regression
    theta=inv(phi_SD'*phi_SD+eye(3)*0.00000001)*phi_SD'*y_SD;

    %Data-Driven errors
    %Generate non-noisy demoinstration using the learned parameters
    [q_l,q_dot_l,u_l] = get_demo(system, time, ts, theta,0,goal);

    % Set the expectation of each risk metrics, the average return of these
    % metrics over large amout of iterations
    Aeloss_avg=	0.002001053;
    L_huber_avg=	1.34644E-05;
    MSLE_avg=	2.34509E-06;
    NMSE_avg=	5.25887E-06;
    RMSE_avg=	0.000122531;
    RMSE_u_avg=	0.015502015;
    NMSE_u_avg=	0.000304375;
    Aerror_avg=	0.006013127;
    MSLE_u_avg=	0.039523525;

    %Calculate the loss between the data from the learned
    %parameters and the nin-noisy demo using the skill parameters
    %MSRE over predicted states
    RMSE = sqrt(((q-q_l))'*(q-q_l))/length(q);

    %MSRE_u over predicted states
    RMSE_u = sqrt(((u-u_l))'*(u-u_l))/length(u);

    %NMSE over predicted action
    NMSE=((q-q_l)'*(q-q_l))/(q'*q);

    %NMSE over predicted action
    NMSE_u=((u-u_l)'*(u-u_l))/(u'*u);

    % %Area between desired path and learned path.

    x = l * cos(q) ;
    y = l * sin(q) ;
    x_l = l * cos(q_l);
    y_l = l * sin(q_l);
    AError= sum(abs((((x-x_l).^2+(y-y_l).^2).^0.5)*ts));

    %Absolute Error Loss over regenirated data
    AEloss=sum(abs(q-q_l))/length(q);

    %Huber Loss
    segma=0.1;
    L=abs(q-q_l);
    for n=1:1:length(L)
        if (L(n) <=  segma)
            L(n)=L(n)^2;
        else
            L(n)=2*segma*(L(n)-segma/2);
        end

    end
    L_huber=sum(L)/length(L);

    %Mean Squared Logarithmic Error (MSLE)
    MSLE=(log(q+ones(length(q),1))-log(q_l+ones(length(q_l),1)))'*(log(q+ones(length(q),1))-log(q_l+ones(length(q_l),1)))/length(q);
    MSLE_u=(log(u+ones(length(u),1))-log(u_l+ones(length(u_l),1)))'*(log(u+ones(length(u),1))-log(u_l+ones(length(u_l),1)))/length(u);

    %Calculate the Learned loss, paramtrized loss
    Leaned_loss = 1.0945*AEloss/Aeloss_avg + 1.7693*L_huber/L_huber_avg + 3*MSLE/MSLE_avg + 1.4336*NMSE/NMSE_avg + 1.4336*RMSE/RMSE_avg + 3*RMSE_u/RMSE_u_avg + 3*NMSE_u/NMSE_u_avg + 1.0945*AError/Aerror_avg + 0.5*MSLE_u/MSLE_u_avg;


    %Calculate the paramter loss
    loss=((sum(abs(theta-theta_opt))/(theta_opt(1)^2+theta_opt(2)^2+theta_opt(3)^2)^0.5))*100;
    % Sotring Data
    storing_list(i,1)=N1;
    storing_list(i,2)=N2;
    storing_list(i,3)=N3;
    storing_list(i,4)=AEloss;
    storing_list(i,5)=L_huber;
    storing_list(i,6)=MSLE;
    storing_list(i,7)=NMSE;
    storing_list(i,8)=RMSE;
    storing_list(i,9)=RMSE_u;
    storing_list(i,10)=NMSE_u;
    storing_list(i,11)=AError;
    storing_list(i,12)=MSLE_u;
    storing_list(i,13)=Leaned_loss;
    storing_list(i,14)=loss;



end


%Demo
    function [q_with_noise,q_with_noise_dot,u_noise,integral_e] = get_demo(system, time, ts, theta,noize,goal)
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
        integral_e=zeros(num_steps + 1, 1);
        integral_e(1)=goal*ts;
        q_with_noise(1)=qi;
        e(1)=goal-qi;
        q_with_noise_dot(1)=q_doti;
        u_noise(1)= theta'*[e(1) integral_e(1) -q_with_noise_dot(1)]';

        for sys_iteration=1:num_steps

            e(sys_iteration)=goal-q_with_noise(sys_iteration);
            integral_e(sys_iteration)=sum(e)*ts;
            phi=[e(sys_iteration) integral_e(sys_iteration) -q_with_noise_dot(sys_iteration)]';
            u_noise(sys_iteration)=(theta'*phi)*(1+noize*randn(1));
            q_with_noise_dot(sys_iteration+1)=q_with_noise_dot(sys_iteration)+((u_noise(sys_iteration)/m*l.^2)-g*cos(q_with_noise(sys_iteration))/l)*ts;
            q_with_noise(sys_iteration+1)=q_with_noise(sys_iteration)+q_with_noise_dot(sys_iteration)*ts;
            q_with_noise_dot(sys_iteration)=q_with_noise_dot(sys_iteration)*(1-noize*rand(1));
            q_with_noise(sys_iteration)=q_with_noise(sys_iteration)*(1-noize*rand(1));
        end
    end
end


