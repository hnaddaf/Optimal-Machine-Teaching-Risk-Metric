# Optimal-Machine-Teaching-Risk-Metric
In this project there is three experiments: (1) for single-link robot arm (2) for double-link robot arm (3) Validation Experiment.

# Machine-Teaching-of-a-single-link-robot

# Overview:
This project investigates the effect of different data-driven loss functions on the machine teaching problem of teaching a single-robot arm a specific skill. It takes the skill parameters and the system variables as input and returns the mean of loss by each risk metric and the mean and standard deviation of the parameter loss associated with all the risk metrics. The process starts by setting up the system and skill parameters. A zero-matrices called General Loss, and Parameter Loss with dimensions of Trials x 11 is initialized to store the results of each trial. Each column of the 11 columns contains the risk loss, and parameter loss associated with each one of the risk metrics.

# How to use:
Run Experiment1.m
In the first of the code the skill parameters for diffrent skills are set as the following:

%Set the skill parameters 

kps=[ 178; 400; 94; 400; 170; 81; 40];

kds=[ 24; 50; 22; 50; 20; 12; 12];

kis=[ 20; 4000; 20; 800; 220; 320; 120];

You can change these PID gains to any policy you prefer. The code then measures the effect of each risk metric and returns the results in an Excel sheets:

1_ experiment1_PARAMETER_LOSS_new.xlsx:

This sheet contains information on the parameter loss for each data-driven risk metric. Parameter loss is the difference between the defined skill parameter and the parameters learned from the data. The data has been selected to minimize the error for each risk metric.

Key details:

Mean and Standard Deviation: These values are provided for the parameter loss.
Trials: The calculations are based on 1000 trials.

1_ experminet1_RISK_LOSS_new.xlsx:
This sheet contains the mean of minimum loss obtained using each one of the data-driven risk metric. 

Trials: The calculations are based on 1000 trials.

# Requirements:
MATLAB


# Machine-Teaching-of-a-Double-link-robot

# Overview:
This project investigates the effect of different data-driven loss functions on the machine teaching problem of teaching a Double-robot arm a specific skill. It takes the skill parameters and the system variables as input and returns the mean of loss by each risk metric and the mean and standard deviation of the parameter loss associated with all the risk metrics. The process starts by setting up the system and skill parameters. A zero-matrices called General Loss, and Parameter Loss with dimensions of Trials x 11 is initialized to store the results of each trial. Each column of the 11 columns contains the risk loss, and parameter loss associated with each one of the risk metrics.

# How to use:
The main code to run is Experiment2.m
In the first of the code, the skill, which is defined as different goals of the end-effector, is set as the following:

% Set the goal joint angles, the target position

goals1=[pi/6; 2*pi/6; 3*pi/6; 4*pi/6; 5*pi/6; pi; 7*pi/6;8*pi/6;9*pi/6;pi/6;pi/4;pi/3;pi/2;2*pi/3;5*pi/6;pi;pi/12;pi/8;pi/6;pi/4;2*pi/6;5*pi/12;pi/2];

goals2=[0;1;-1;0;0.5;-0.6;1.2;pi;pi/2;pi/12;pi/8;pi/6;pi/4;2*pi/6;5*pi/12;pi/2;pi/6;pi/4;pi/3;pi/2;2*pi/3;5*pi/6;pi];


You can change these goals to any targets you prefer. The code then measures the effect of each risk metric and returns the results in an Excel sheets:

1_ experminet1_PARAMETER_LOSS_new.xlsx:

This sheet contains information on the parameter loss for each data-driven risk metric. Parameter loss is the difference between the defined skill parameter and the parameters learned from the data. The data has been selected to minimize the error for each risk metric.

Key details:

Mean and Standard Deviation: These values are provided for the parameter loss.
Trials: The calculations are based on 1000 trials.

1_ experminet1_RISK_LOSS_new.xlsx:
This sheet contains the mean of minimum loss obtained using each one of the data-driven risk metric. 

Trials: The calculations are based on 1000 trials.

# Requirements:
MATLAB

# Validation Experiment

# Overview:
This part validate the methodlodgy by applying it to the optemal data-driven machine teaching for PD controller policy. 

# How to use:
Run ValidationExperiment.m
The code then measures parameter loss and optimal data-driven loss and returns the results in an Excel sheets:

# Requirements:
MATLAB
