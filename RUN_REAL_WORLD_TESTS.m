%% RUNS EKF on Collected Data
close all; clear;
addpath(genpath('experiments'))
addpath(genpath('Segmented_data'))
addpath(genpath('EKF'))

load('SET1')
%load('SET2') % has some problems with getting correct bias 
%load('SET3')

load('Mag_Calibration_Data.mat')

%% Settings
% Get Settings struct
EKF_settings = get_settings_experimentsD1();
% Get initial state from filter
X_initial = [[0 0 0]... % Initial Position
                 Data_N.EKF.Vel(:,1)' ... % Initial Velocity
                 dcm2euler(quat2rotm(Data_N.EKF.ORIq(:,1)')) ... % Initial Attitude (euler)
                 -Data_N.EKF.ACCB(:,1)' ... % Initial Acc Bias
                 -Data_N.EKF.GYROB(:,1)' ]';  % Initial Gyro Bias
I = eye(3); 
P_initial = blkdiag(Data_N.EKF.POSunc(:,1).^2.*I, ...
            Data_N.EKF.Velunc(:,1).^2.*I, ...
            Data_N.EKF.ORIunc(:,1).^2.*I, ...
            Data_N.EKF.ACCBunc(:,1).^2.*I, ...
            Data_N.EKF.GYROBunc(:,1).^2.*I); 

%%
% Set Gravity Vector from device
EKF_settings.g = -[0 0 Data_N.EKF.GRAVITY(end)];
% Set magnetic Field vector from device. This should likley be verify with a table         
EKF_settings.MagneticField = Data_N.EKF.EARTH_MAGNETIC_VEC(:,1)';
% Calibrate Magnetometer settings
Data_N.IMU.mag = ((Data_N.IMU.mag - Mag_Calibration.HardIron')'*Mag_Calibration.SoftIron)';

% Rename data variable
Data = Data_N;
% Add pressurments under correct name for EKF
Data.PRESS_SEN = Data_N.pressure;
% Remove the atomshpheric pressure, this varies a bit but it is assumed
% constant
Data.PRESS_SEN.press = (Data.PRESS_SEN.press - 1026);
% Calculate Zero Velocity Update
EKF_settings.gamma = 300;
EKF_settings.window_size = 50; 
[T_zupt, zupt] = ZUPT_CAL([Data_N.IMU.acc; Data_N.IMU.gyro],EKF_settings,'GLRT'); 
% Add Zupt to Data File
Data.ZUPT.t = Data.IMU.t;
Data.ZUPT.ZUPT = zupt;

% Add Orientation estiamte Data from INS sesnor
Data.COMPASS.t = Data_N.EKF.t;
Data.COMPASS.Data = flip(quat2eul(Data_N.EKF.ORIq')');
% Orientation estimate from compass, sliglty different to INS sensor
% Corrently the INS sesnor values are used
% The magnetic vector could be an error source
% Using the INS sesnors attitude measurment generally perfomed better
Data.COMPASS2.data = [];
Data.COMPASS2 = get_compass(Data.IMU.acc,Data.IMU.mag,Data.IMU.t,11,EKF_settings);

% Remove bad GPS data
% GPS data is given even if no Fix is present
% Threshold is used to remove data not with the specific bound
thresshold = 5;
Data.GNSS.t = Data.GNSS.t((max(Data.GNSS.accuracy) < thresshold));
Data.GNSS.pos = Data.GNSS.pos(:,(max(Data.GNSS.accuracy) < thresshold));
Data.GNSS.accuracy = Data.GNSS.accuracy(:,(max(Data.GNSS.accuracy) < thresshold));
%Data.GNSS.accuracy = sqrt(Data.GNSS.accuracy);
% Adjust down axis, sicne sea lvl is not neccesarily at zero point
X_initial(3) =  - mean(Data.GNSS.pos(3,:)) + 0.3; 
Data.GNSS.pos(3,:)= Data.GNSS.pos(3,:) - mean(Data.GNSS.pos(3,:)) + 0.3;
% Find position where it is dropped
[~,n] = max(diff(Data.GNSS.t));
time_drop = Data.GNSS.t(n);
time_pick = Data.GNSS.t(n+1);
%% Set Settings
EKF_settings.sigma_acc = 0.60; % [m/s^2] %0.25 %0.30
EKF_settings.sigma_gyro =0.15*pi/180; % [rad/s] %=0.05 &0.048*pi/180; %0.048
EKF_settings.sigma_mag = 0.05; % Measurment noise Mag
EKF_settings.sigma_zero_velocity = 0.000000001; % ZUPT Noise
settings.acc_bias_driving_noise=0.000001; % Acc bias change rate
settings.gyro_bias_driving_noise=0.00001*pi/180; % Gyro bias change rate
EKF_settings.compass_cutoff = 0.58; %error rejection for compass
EKF_settings.compass_pressure_cutoff = 50; %Stop reading since boat will effect meas
%% Prep variables
X = X_initial;
P = P_initial;
X_n = X;
P_n = P;
DataS = Data;
DataS.IMU.t = DataS.IMU.t(1:floor(length(Data.IMU.t)/3*2));

biases = X(10:15);
biases_p = diag(P(10:15,10:15));
%% Perform Iterative model
for o=1:10 %1:15
    
    % Kalman filtering
    outdata= EKF(X_n, P_n,DataS,EKF_settings);
    
    % Smoothing
    soutdata = RTS2(outdata,EKF_settings);
    
    % Update Biases
    X_n = soutdata.X(:,1); 
    P_n = (((soutdata.P(:,:,1))));
    X_n = [X(1:9);  X_n(10:end,1)]; 
    P_n = blkdiag((P(1:9,1:9)), ((P_n(10:end,10:end))));
    
    Position_itt = [outdata.X(1:3,end)' o]
    figure(13);
    plot(outdata.X(10:12,:)')
    title('Acc Biases')
    figure(14);
    plot(outdata.X(13:15,:)')
    title('Gyro Biases')
    
    biases(:,o+1) = X_n(10:15)';
    biases_p(:,o+1) = diag(P_n(10:15,10:15)');
    
end
%% Get results

outdata = EKF(X_n, P_n,Data,EKF_settings);
outdata_n = EKF(X, P,Data,EKF_settings);

pos_estimate_n = outdata_n.X(1:3,floor(length(Data.IMU.t)/3*2))';
pos_estimate_conf_n = sqrt(outdata_n.diag_P(1:3,floor(length(Data.IMU.t)/3*2)))';
notsmoothed_position = [pos_estimate_n pos_estimate_conf_n];

pos_estimate = outdata.X(1:3,floor(length(Data.IMU.t)/3*2))';
pos_estimate_conf = sqrt(outdata.diag_P(1:3,floor(length(Data.IMU.t)/3*2)))';
smoothed_position = [pos_estimate pos_estimate_conf];

% Drop/pick up position
drop = [Data.GNSS.pos(:,n)' Data.GNSS.accuracy(1,n) Data.GNSS.accuracy(:,n)'];
pick = [Data.GNSS.pos(:,n+1)' Data.GNSS.accuracy(1,n+1) Data.GNSS.accuracy(:,n+1)'];

RESULTS = [pos_estimate_n pos_estimate_conf_n;pos_estimate pos_estimate_conf; drop; pick];

TYPE = ["EKF" "Iterative EKF" "Drop-off Position" "Pick up position"]';
NORTH_POSITION = RESULTS(:,1);
EAST_POSITION = RESULTS(:,2);
DOWN_POSITION = RESULTS(:,3);
NORTH_UNC_POSITION = RESULTS(:,4);
EAST_UNC_POSITION = RESULTS(:,5);
DOWN_UNC_POSITION = RESULTS(:,6);

T = table(TYPE,NORTH_POSITION,EAST_POSITION,DOWN_POSITION,NORTH_UNC_POSITION,EAST_UNC_POSITION,DOWN_UNC_POSITION)


plot_real_trajectory(soutdata,time_drop-50,Data, drop, pick,'Smoothed Path Iterative Version', 1)

%% Path Estimate
figure(10);
sgtitle('Position Result Dataset 1 - Iterative Bias Update')
% Plot
subplot(3,1,1)
plot(outdata.t,outdata.X(1,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(1,:),'.k')

title('North')
xlabel('time (s)')
ylabel('pos (m)')
subplot(3,1,2)
plot(outdata.t,outdata.X(2,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(2,:),'.k')

title('East')
xlabel('time (s)')
ylabel('pos (m)')
subplot(3,1,3)
plot(outdata.t,outdata.X(3,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(3,:),'.k')
title('Down')
xlabel('time (s)')
ylabel('pos (m)')

figure(11);
sgtitle('Position Result Dataset 1 - EKF')
% Plot
subplot(3,1,1)
plot(outdata_n.t,outdata_n.X(1,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(1,:),'.k')

title('North')
xlabel('time (s)')
ylabel('pos (m)')
subplot(3,1,2)
plot(outdata_n.t,outdata_n.X(2,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(2,:),'.k')

title('East')
xlabel('time (s)')
ylabel('pos (m)')
subplot(3,1,3)
plot(outdata_n.t,outdata_n.X(3,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(3,:),'.k')
title('Down')
xlabel('time (s)')
ylabel('pos (m)')

 
%% Biases Improvment
figure(12);
sgtitle('Initial Accelerometer Bias Estimate - Set 1')
% Plot
subplot(3,1,1)
plot(biases(1,:))
hold on 
plot(biases(1,:) + 3*sqrt(biases_p(1,:)), 'r');
plot(biases(1,:) - 3*sqrt(biases_p(1,:)), 'r');
hold off
title('North')
xlabel('Iteration')
ylabel('Bias (m/s^2)')
subplot(3,1,2)
plot(biases(2,:))
hold on 
plot(biases(2,:) + 3*sqrt(biases_p(2,:)), 'r');
plot(biases(2,:) - 3*sqrt(biases_p(2,:)), 'r');
hold off  
title('East')
xlabel('Iteration')
ylabel('Bias estimate (m/s^2)')
subplot(3,1,3)
plot(biases(3,:))
hold on 
plot(biases(3,:) + 3*sqrt(biases_p(3,:)), 'r');
plot(biases(3,:) - 3*sqrt(biases_p(3,:)), 'r');
hold off
title('Down')
xlabel('Iteration')
ylabel('Bias (m/s^2)')

figure(13);
sgtitle('Initial Gyroscope Bias Estimate - Set 1')
% Plot
subplot(3,1,1)
plot(biases(4,:))
hold on 
plot(biases(4,:) + 3*sqrt(biases_p(4,:)), 'r');
plot(biases(4,:) - 3*sqrt(biases_p(4,:)), 'r');
hold off
title('North')
xlabel('Iteration')
ylabel('Bias (rad/s)')
subplot(3,1,2)
plot(biases(5,:))
hold on 
plot(biases(5,:) + 3*sqrt(biases_p(5,:)), 'r');
plot(biases(5,:) - 3*sqrt(biases_p(5,:)), 'r');
hold off 
title('East')
xlabel('Iteration')
ylabel('Bias estimate (rad/s)')
subplot(3,1,3)
plot(biases(6,:))
hold on 
plot(biases(6,:) + 3*sqrt(biases_p(6,:)), 'r');
plot(biases(6,:) - 3*sqrt(biases_p(6,:)), 'r');
hold off
title('Down')
xlabel('Iteration')
ylabel('Bias (rad/s)')

%% Compass vs INS
% This relies highly on the magnetic vector 

figure(15);
sgtitle('Attitude from own compass VS Attitude INS sensors')
% Plot
subplot(3,1,1)
h(1) = plot(Data.COMPASS.t,Data.COMPASS.Data(1,:),'r');
hold on 
h(2) = plot(Data.COMPASS2.t,Data.COMPASS2.Data(1,:),'b');
hold off
legend(h,{'INS COMPASS','MAG COMPASS'})
title('Roll')
xlabel('t(s)')
ylabel('rad')
axis([0 Data.COMPASS.t(end) -pi pi])
subplot(3,1,2)
h(1) = plot(Data.COMPASS.t,Data.COMPASS.Data(2,:),'r');
hold on 
h(2) = plot(Data.COMPASS2.t,Data.COMPASS2.Data(2,:),'b');
hold off
legend(h,{'INS COMPASS','MAG COMPASS'})
title('Pitch')
xlabel('t(s)')
ylabel('rad')
axis([0 Data.COMPASS.t(end) -pi pi])
subplot(3,1,3)
h(1) = plot(Data.COMPASS.t,Data.COMPASS.Data(3,:),'r');
hold on 
h(2) = plot(Data.COMPASS2.t,Data.COMPASS2.Data(3,:),'b');
hold off
legend(h,{'INS COMPASS','MAG COMPASS'})
title('Yaw')
xlabel('t(s)')
ylabel('rad')
axis([0 Data.COMPASS.t(end) -pi pi])
