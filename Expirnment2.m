% Set the goal joint angles, the target position
goals1=[pi/6; 2*pi/6; 3*pi/6; 4*pi/6; 5*pi/6; pi; 7*pi/6;8*pi/6;9*pi/6;pi/6;pi/4;pi/3;pi/2;2*pi/3;5*pi/6;pi;pi/12;pi/8;pi/6;pi/4;2*pi/6;5*pi/12;pi/2];
goals2=[0;1;-1;0;0.5;-0.6;1.2;pi;pi/2;pi/12;pi/8;pi/6;pi/4;2*pi/6;5*pi/12;pi/2;pi/6;pi/4;pi/3;pi/2;2*pi/3;5*pi/6;pi];
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

for systems=1:length(goals1)% Repteate the experiment of the diffrent skills
    %Set the initial values of the paramters for the Paramtrized loss
    paramters=[1 1 1 1 1 1 1 1 1];
    %Set the time parameters
    sampling_time = 0.005;
    time = 2;
    ts = sampling_time;
    num_steps=time/sampling_time;
    % Defined the robot geometry
    m1 = 1;
    l1 = 1;
    m2 = 1;
    l2 = 0.5;
    g = 9.81;

    
    system=[m1,l1,m2,l2,g];

    % Set PID gains
    kp1 = 350;
    kp2 = 200;
    ki1 = 300;
    ki2 = 200;
    kd1 = 60;
    kd2 = 10;
    %Set skill parameters
    theta1_star=[kp1; ki1; kd1];
    theta2_star=[kp2; ki2; kd2];
    %Set the target joint angles
    q1f=goals1(systems);
    q2f=goals2(systems);


    trials=1000;
    OD1=zeros(trials,11);% Set the Paramter Loss matrics
    OD2=zeros(trials,11);% Set the General Loss matrics
    for n=1:trials % Rpeate the experiment for 1000 trials

        % Generate a noisy demonistration with noise with the mean of 5% of
        % the states and actions value
        [q1_noise q2_noise q1dot_noise q2dot_noise int_e1_noise int_e2_noise u1_noise u2_noise] = getdemo_with_noise(time,ts,theta1_star,theta2_star,system,q1f,q2f,0.05);
        % Generate non noisy demonistration
        [q1 q2 q1dot q2dot int_e1 int_e2 u1 u2] = getdemo_with_noise(time,ts,theta1_star,theta2_star,system,q1f,q2f,0.0);

        
        iteration=1000;
        storing_list=zeros(iteration,17);% Create the sorting list as zero matrix

        for i = 1:iteration % Get 1000 line of data into sorting list

            %Data Selection, generate six data points in the range of
            %number of samples
            while true
                N11=round((time/ts)*rand(1));
                N21=round((time/ts)*rand(1));
                N31=round((time/ts)*rand(1));
                if (N11~=0)&&(N21 ~= 0)&&(N31~=0)
                    break
                end
            end
            while true
                N12=round((time/ts)*rand(1));
                N22=round((time/ts)*rand(1));
                N32=round((time/ts)*rand(1));
                if (N12~=0)&&(N22 ~= 0)&&(N32~=0)
                    break
                end
            end
            % Use the random data to select the data out of the noisy demo
            % create the feature vectors
            y_SD1=[u1_noise(N11);u1_noise(N21);u1_noise(N31)];
            phi_SD1=[q1f-q1_noise(N11) q1f-q1_noise(N21) q1f-q1_noise(N31);int_e1_noise(N11) int_e1_noise(N21) int_e1_noise(N31) ;-q1dot_noise(N11) -q1dot_noise(N21)  -q1dot_noise(N31)]';
            y_SD2=[u2_noise(N12);u2_noise(N22);u2_noise(N32)];
            phi_SD2=[q2f-q2_noise(N12) q2f-q2_noise(N22) q2f-q2_noise(N32);int_e2_noise(N12) int_e2_noise(N22) int_e2_noise(N32) ;-q2dot_noise(N12) -q2dot_noise(N22)  -q2dot_noise(N32)]';

            %Learning, Ridge regression
            theta1=inv(phi_SD1'*phi_SD1+eye(3)*0.00000001)*phi_SD1'*y_SD1;
            theta2=inv(phi_SD2'*phi_SD2+eye(3)*0.00000001)*phi_SD2'*y_SD2;

            %Data-Driven errors
            %Generate non-noisy demoinstration using the learned parameters
            [q1_l q2_l q1dot_l q2dot_l int_e1_l int_e2_l u1_l u2_l] = getdemo_with_noise(time,ts,theta1,theta2,system,q1f,q2f,0.0);
            

            %Calculate the loss between the data from the learned
            %parameters and the nin-noisy demo using the skill parameters

            %Absolute Error Loss over regenirated data
            AEloss=(sum(abs(q1-q1_l))/length(q1))+(sum(abs(q2-q2_l))/length(q2));

            %Huber Loss
            segma=0.01;
            L1=abs(q1-q1_l);
            for ni=1:1:length(L1)
                if (L1(ni) <=  segma)
                    L1(ni)=L1(ni)^2;
                else
                    L1(ni)=2*segma*(L1(ni)-segma/2);
                end

            end

            L2=abs(q2-q2_l);
            for ni=1:1:length(L2)
                if (L2(ni) <=  segma)
                    L2(ni)=L2(ni)^2;
                else
                    L2(ni)=2*segma*(L2(ni)-segma/2);
                end

            end
            L_huber=(sum(L1)+sum(L2))/length(L1);


            %Mean Squared Logarithmic Error (MSLE) over leared path
            MSLE1=(log(q1+ones(length(q1),1))-log(q1_l+ones(length(q1_l),1)))'*(log(q1+ones(length(q1),1))-log(q1_l+ones(length(q1_l),1)))/length(q1);
            MSLE2=(log(q2+ones(length(q2),1))-log(q2_l+ones(length(q2_l),1)))'*(log(q2+ones(length(q2),1))-log(q2_l+ones(length(q2_l),1)))/length(q2);
            MSLE=MSLE1+MSLE2;

            %NMSE over predicted states
            NMSE=(((q1-q1_l)'*(q1-q1_l))/(q1'*q1))+((q2-q2_l)'*(q2-q2_l))/(q2'*q2);

            %RMSE over predicted states
            RMSE = (sqrt(((q1-q1_l))'*(q1-q1_l))/length(q1))+sqrt(((q2-q2_l))'*(q2-q2_l))/length(q2);

            %RMSE over predicted states
            RMSE_u = (sqrt(((u1-u1_l))'*(u1-u1_l))/length(u1))+(sqrt(((u2-u2_l))'*(u2-u2_l))/length(u2));

            %NMSE over predicted action
            NMSE_u=(((u1-u1_l)'*(u1-u1_l))/(u1'*u1))+(((u2-u2_l)'*(u2-u2_l))/(u2'*u2));

            % %Area between desired path and learned path.
            x = l1 * cos(q1) + l2 * cos(q1 + q2);
            y = l1 * sin(q1) + l2 * sin(q1 + q2);
            x_l = l1 * cos(q1_l) + l2 * cos(q1_l + q2_l);
            y_l = l1 * sin(q1_l) + l2 * sin(q1_l + q2_l);

            AError= sum(abs((((x-x_l).^2+(y-y_l).^2).^0.5)*ts));


            %Mean Squared Logarithmic Error (MSLE) over leared action
            MSLE1_u=(log(u1+ones(length(u1),1))-log(u1_l+ones(length(u1_l),1)))'*(log(u1+ones(length(u1),1))-log(u1_l+ones(length(u1_l),1)))/length(u1);
            MSLE2_u=(log(u2+ones(length(u2),1))-log(u2_l+ones(length(u2_l),1)))'*(log(u2+ones(length(u2),1))-log(u2_l+ones(length(u2_l),1)))/length(u2);
            MSLE_u=MSLE1_u+MSLE2_u;
         
            %Calculate the Learned loss, paramtrized loss
            Leaned_loss =     2.0751*AEloss/Aeloss_avg + 2.2797*L_huber/L_huber_avg + 2.2932*MSLE/MSLE_avg + 2.2508*NMSE/NMSE_avg + 2.3016*RMSE/RMSE_avg + 0.4615*RMSE_u/RMSE_u_avg + 0.5935*NMSE_u/NMSE_u_avg ;
    

            %Calculate the paramter loss
            loss1=((sum(abs(theta1_star-theta1))/(theta1_star(1)^2+theta1_star(2)^2+theta1_star(3)^2)^0.5))*100;
            loss2=((sum(abs(theta2_star-theta2))/(theta2_star(1)^2+theta2_star(2)^2+theta2_star(3)^2)^0.5))*100;
            loss=(loss1+loss2)/2;


            % Sotring Data
            storing_list(i,1)=N11;
            storing_list(i,2)=N21;
            storing_list(i,3)=N31;
            storing_list(i,4)=N12;
            storing_list(i,5)=N22;
            storing_list(i,6)=N32;
            storing_list(i,7)=AEloss;
            storing_list(i,8)=L_huber;
            storing_list(i,9)=MSLE;
            storing_list(i,10)=NMSE;
            storing_list(i,11)=RMSE;
            storing_list(i,12)=RMSE_u;
            storing_list(i,13)=NMSE_u;
            storing_list(i,14)=AError;
            storing_list(i,15)=MSLE_u;
            storing_list(i,16)=Leaned_loss;
            storing_list(i,17)=loss;

        end


        % Locate the minimum loss optained by each risk metrics and return its value and location 
        [loss1 ,location_of_Data1]= min(storing_list(:,7));
        [loss2 ,location_of_Data2]= min(storing_list(:,8));
        [loss3 ,location_of_Data3]= min(storing_list(:,9));
        [loss4 ,location_of_Data4]= min(storing_list(:,10));
        [loss5 ,location_of_Data5]= min(storing_list(:,11));
        [loss6 ,location_of_Data6]= min(storing_list(:,12));
        [loss7 ,location_of_Data7]= min(storing_list(:,13));
        [loss8 ,location_of_Data8]= min(storing_list(:,14));
        [loss9 ,location_of_Data9]= min(storing_list(:,15));
        [loss10 ,location_of_Data10]= min(storing_list(:,16));
        [loss11 ,location_of_Data11]= min(storing_list(:,17));

        % Save the parameter loss at the location of minimum risk loss in L 
        L1=storing_list(location_of_Data1,17);
        L2=storing_list(location_of_Data2,17);
        L3=storing_list(location_of_Data3,17);
        L4=storing_list(location_of_Data4,17);
        L5=storing_list(location_of_Data5,17);
        L6=storing_list(location_of_Data6,17);
        L7=storing_list(location_of_Data7,17);
        L8=storing_list(location_of_Data8,17);
        L9=storing_list(location_of_Data9,17);
        L10=storing_list(location_of_Data10,17);
        L11=min(storing_list(:,17));

        % Learn the paramters for the Paramtrized loss
        L=[L1;L2;L3;L4;L5;L6;L7;L8;L9];

        mean= sum(L)/9;
        for in=1:9

            paramters(in)=paramters(in)-0.0005*(L(in)-mean);

        end

        % append the parameter loss at the location of minimum risk loss to
        % the Parameter loss matrix
        OD1(n,1)=L1;
        OD1(n,2)=L2;
        OD1(n,3)=L3;
        OD1(n,4)=L4;
        OD1(n,5)=L5;
        OD1(n,6)=L6;
        OD1(n,7)=L7;
        OD1(n,8)=L8;
        OD1(n,9)=L9;
        OD1(n,10)=L10;
        OD1(n,11)=L11;

        % append the minimum risk loss to the General loss matrix
        OD2(n,1)=loss1;
        OD2(n,2)=loss2;
        OD2(n,3)=loss3;
        OD2(n,4)=loss4;
        OD2(n,5)=loss5;
        OD2(n,6)=loss6;
        OD2(n,7)=loss7;
        OD2(n,8)=loss8;
        OD2(n,9)=loss9;
        OD2(n,10)=loss10;
        OD2(n,11)=loss11;


    end

    % Create a list to store the mean and the standard deviation of the
    % parameter loss assosiated with each risk metric
    list2={'Loss Metrics', 'mean of parameter loss','Standard Deviation';
        'AEloss',sum(OD1(:,1))/trials ,std(OD1(:,1));
        'L_huber',sum(OD1(:,2))/trials ,std(OD1(:,2)) ;
        'MSLE',sum(OD1(:,3))/trials ,std(OD1(:,3));
        'NMSE',sum(OD1(:,4))/trials ,std(OD1(:,4));
        'RMSE',sum(OD1(:,5))/trials ,std(OD1(:,5));
        'RMSE_u',sum(OD1(:,6))/trials ,std(OD1(:,6)) ;
        'NMSE_u',sum(OD1(:,7))/trials ,std(OD1(:,7)) ;
        'AError',sum(OD1(:,8))/trials ,std(OD1(:,8));
        'MSLE_u',sum(OD1(:,9))/trials ,std(OD1(:,9));
        'Leaned_loss',sum(OD1(:,10))/trials ,std(OD1(:,10));
        'Parameter loss',sum(OD1(:,11))/trials ,std(OD1(:,11))}
     % Create a list to store the mean of the risk loss assosiated with each risk metric
     list3={'Loss Metrics', 'mean of risk loss';
        'AEloss',sum(OD2(:,1))/trials;
        'L_huber',sum(OD2(:,2))/trials ;
        'MSLE',sum(OD2(:,3))/trials ;
        'NMSE',sum(OD2(:,4))/trials ;
        'RMSE',sum(OD2(:,5))/trials ;
        'RMSE_u',sum(OD2(:,6))/trials;
        'NMSE_u',sum(OD2(:,7))/trials ;
        'AError',sum(OD2(:,8))/trials;
        'MSLE_u',sum(OD2(:,9))/trials;
        'Leaned_loss',sum(OD2(:,10))/trials ;
        'Parameter loss',sum(OD2(:,11))/trials}

     % Save the lists to an excel file
    filename = 'experminet2_Parameter_Loss.xlsx';
    writecell(list2,filename,'Sheet',systems,'Range','D1')

    filename = 'experminet2_Risk_Loss.xlsx';
    writecell(list3,filename,'Sheet',systems,'Range','D1')

end