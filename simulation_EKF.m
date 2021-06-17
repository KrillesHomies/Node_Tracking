%% This script runs the the Kalman filter on the simulated version
% All complemtary functions are in the simulation folder,
% the data can be found there too. New data sets can also be easily created
% Last edited 2021-06-17 by Kristoffer Lindve
%%
close all; clear all;
addpath(genpath('EKF'))
addpath(genpath('simulation'))

% load data - Data set is choose by uncommenting it
% load('Data_set20.mat')
% load('Data_set30.mat')
load('Data_set40.mat')

% Get settings
EKF_settings = get_sim_noise_settings();

%% Inital variables
I = eye(3); 
X_initial = [0 0 0 ... % Initial Position
                 0 0 0 ... % Initial Velocity
                 0 0 0 ... % Initial Attitude (euler)
                 0 0 0 ... % Initial Acc Bias
                 0 0 0]';  % Initial Gyro Bias 
P_initial = blkdiag(EKF_settings.factp(1)^2*I, ... % Initial Uncertainity Position
            EKF_settings.factp(2)^2*I, ...% Initial Uncertainity Velocity
            diag(EKF_settings.factp(3:5)).^2, ... % Initial Uncertainity Attitude
            EKF_settings.factp(6)^2*I, ... % Initial Uncertainity Acc Bias
            EKF_settings.factp(7)^2*I);   % Initial Uncertainity Gyro Bias
X = X_initial; 
P = P_initial;
%%
%Create Randomized color scheme
colors = linspace(0.1,1.0,size(Dataset,2)).*[1; 1; 1];
colors(1:prod(size(colors))) = colors(randperm(prod(size(colors))));


half_way_point = floor(length(Dataset{1}.Data_N.IMU.t)/2);
half_way_point_gnss = floor(length(Dataset{1}.Data_NN.GNSS.t(:))/2);

Data_vector_EKF = [];
Data_vector_itt = [];

% Used for stats
stats_ekf = [];
stats_itt = [];

for j=1:size(Dataset,2)
    
    Data = Dataset{j}.Data_N;
    % Remove GPS points that under the water surface 
    i = (Dataset{j}.Data_NN.GNSS.pos(3,:) < 0);
    Data.GNSS.pos = Data.GNSS.pos(:,i);
    Data.GNSS.t = Data.GNSS.t(i);
    % Get ZUPT values and add them to the data array
    EKF_settings = get_sim_noise_settings();
    EKF_settings.gamma = 5;
    EKF_settings.window_size = 10;
    [T_zupt, zupt] = ZUPT_CAL([Data.IMU.acc; Data.IMU.gyro],EKF_settings,'GLRT');
    Data.ZUPT.t = Data.IMU.t;
    Data.ZUPT.ZUPT = zupt;    
    % Get compass measurments
    Data.COMPASS.t = [];
    Data.COMPASS.data = [];
    Data.COMPASS = get_compass(Data.IMU.acc,Data.IMU.mag,Data.IMU.t,11,EKF_settings);
    
    % Set settings
    EKF_settings.sigma_acc = 0.50 ; % 60 0.20 
    EKF_settings.sigma_gyro = 0.20*pi/180 ;  % 34 0.10*pi/180
    EKF_settings.sigma_mag = 0.25;  
    EKF_settings.compass_cutoff = 0.0001;
    EKF_settings.compass_pressure_cutoff = 50;
    
    
    % EKF on entire dataset
    Data_vector_EKF{j} = EKF(X, P,Data,EKF_settings);
    
    % Iterativly improve bias
    X_n = X;
    P_n = P;
    DataS = Data;
    % removes half the timestamp which will cause the filter to ignore half the data 
    DataS.IMU.t = DataS.IMU.t(1:half_way_point); 
    for o=1:10
        % EKF on half the data
        outdata2 = EKF(X_n, P_n,DataS,EKF_settings);
        % Smoothing
        soutdata = RTS2(outdata2,EKF_settings);
        % Extract only bias
        X_n = [X(1:9);  soutdata.X(10:end,1)]; 
        P_n = blkdiag((P(1:9,1:9)), soutdata.P(10:end,10:end,1));
%         X_n = (soutdata.X(:,1)+X_n)/(2); 
%         P_n = ((((((soutdata.P(:,:,1))))) + (P_n))/(2));
%         X_n = [X(1:9);  X_n(10:end,1)]; 
%         P_n = blkdiag((P(1:9,1:9)), ((P_n(10:end,10:end))));
    end
    Data_vector_itt{j} = EKF(X_n, P_n,Data,EKF_settings);
    
    plot_trajectory(Data_vector_EKF{j}, Dataset{j}.Data_NN,'EKF',10,colors(:,j))
    plot_trajectory(Data_vector_itt{j}, Dataset{j}.Data_NN,'Iterative EKF',11,colors(:,j))
    
    plot_pos_with_bounds(Data_vector_EKF{j}, Dataset{j},half_way_point,'Position Estimate: EKF', 12,colors(:,j))
    plot_pos_with_bounds(Data_vector_itt{j}, Dataset{j},half_way_point,'Position Estimate: Iterative EKF', 13,colors(:,j))
    
    stats_ekf(j) = norm(abs(Dataset{j}.X_r(1:3,half_way_point) - Data_vector_EKF{j}.X(1:3,half_way_point)));
    stats_itt(j) = norm(abs(Dataset{j}.X_r(1:3,half_way_point) - Data_vector_itt{j}.X(1:3,half_way_point)));
end

version = ["EKF"; "Iterative EKF"];
average_error = [mean(stats_ekf); mean(stats_itt)];
std_of_error = [std(stats_ekf); std(stats_itt)];
max_error = [max(stats_ekf); max(stats_itt)];
min_error = [min(stats_ekf); min(stats_itt)];
Distance_NE = [norm(Dataset{j}.X_r(1:2,half_way_point)); norm(Dataset{j}.X_r(1:2,half_way_point))];
Distance_NED = [norm(Dataset{j}.X_r(1:3,half_way_point)); norm(Dataset{j}.X_r(1:3,half_way_point))];
T = table(version,average_error,std_of_error,max_error,min_error,Distance_NE,Distance_NED)
