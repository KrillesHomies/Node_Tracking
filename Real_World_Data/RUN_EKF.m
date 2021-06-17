
close all; clear;
addpath Forsok_1
addpath Forsok_1\Processed_Data
addpath Forsok_1\Cut_Data
%load('D1.mat')
load('2021-04-28_40m_4f_long')
load('Mag_Calibration_Data.mat')
% load('D4.mat')
% Get initial state 
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
% Read Settings  
EKF_settings = get_settings_experimentsD1();
EKF_settings.MagneticField = Data_N.EKF.EARTH_MAGNETIC_VEC(:,end)';
EKF_settings.MagneticField = EKF_settings.MagneticField/norm(EKF_settings.MagneticField)*Mag_Calibration.MagneticFieldStrength;
EKF_settings.g = -[0 0 Data_N.EKF.GRAVITY(end)];
% Calibrate Mag values
Data_N.IMU.mag = ((Data_N.IMU.mag - Mag_Calibration.HardIron')'*Mag_Calibration.SoftIron)';
% Trasnfer data
Data = Data_N;
% Fix pressure sensor part
Data.PRESS_SEN = Data_N.pressure;
Data.PRESS_SEN.press = (Data.PRESS_SEN.press - 1026);
% Zupt
EKF_settings.gamma = 70000;
EKF_settings.window_size = 50; %50 
[T_zupt, zupt] = ZUPT_CAL([Data_N.IMU.acc; Data_N.IMU.gyro],EKF_settings,'GLRT'); 
% zupt(:) = 0;
% zupt(9479:end) = 1;
%%
Data.ZUPT.t = Data.IMU.t;
Data.ZUPT.ZUPT = zupt;

Data.COMPASS.t = Data_N.EKF.t;
Data.COMPASS.Data = flip(quat2eul(Data_N.EKF.ORIq')');

% Gravity eximtation
EKF_settings.window_size = 10;
EKF_settings.gamma = 0.01;
[T, MAG_T] = ZUPT_CAL([Data_N.IMU.acc; Data_N.IMU.gyro],EKF_settings,'MAG');
Data.GRAVITY.t = Data.IMU.t;
Data.GRAVITY.MAG = MAG_T;
% Remove bad GPS data
thresshold = 5;
Data.GNSS.t = Data.GNSS.t((max(Data.GNSS.accuracy) < thresshold));
Data.GNSS.pos = Data.GNSS.pos(:,(max(Data.GNSS.accuracy) < thresshold));
Data.GNSS.accuracy = Data.GNSS.accuracy(:,(max(Data.GNSS.accuracy) < thresshold));
%%
Data.COMPASS.t = Data_N.EKF.t;
Data.COMPASS.Data = flip(quat2eul(Data_N.EKF.ORIq')');
Data.COMPASS.t = [];
%Data.COMPASS.data = [];
%Data.COMPASS = get_compass(Data.IMU.acc,Data.IMU.mag,Data.IMU.t,11,EKF_settings);
EKF_settings.sigma_acc = 0.5; % [m/s^2] %0.25 %0.30
EKF_settings.sigma_gyro =0.15*pi/180; % [rad/s] %=0.05 &0.048*pi/180; %0.048
EKF_settings.sigma_mag = 0.05; % Measurment noise Mag
EKF_settings.sigma_zero_velocity = 0.000000001;
settings.acc_bias_driving_noise=0.0000000001; 
settings.gyro_bias_driving_noise=0.0000000001*pi/180; 
EKF_settings.compass_cutoff = 0.58;
EKF_settings.compass_pressure_cutoff = 50;

X = X_initial;
P = P_initial;
X_n = X;
P_n = P;
DataS = Data;
DataS.IMU.t = DataS.IMU.t(1:floor(length(Data.IMU.t)/3*2));
%%
for o=1:5 %1:15
    
    outdata= EKF(X_n, P_n,DataS,EKF_settings);
    soutdata = RTS2(outdata,EKF_settings);

    X_n = (soutdata.X(:,1)+0*X_n)/(1); 
    P_n = (((abs(diag(diag(soutdata.P(:,:,1))))) + 0*(P_n))/(1));
%         X_n = [X(1:6);  X_n(7:end,1)]; 
%         P_n = blkdiag((P(1:6,1:6)), ((P_n(7:end,7:end))));
    X_n = [X(1:9);  X_n(10:end,1)]; 
    P_n = blkdiag((P(1:9,1:9)), ((P_n(10:end,10:end))));
    [outdata.X(1:3,end)' o]
    figure(13);
    plot(outdata.X(10:12,:)')
    title('Acc Biases')
    figure(14);
    plot(outdata.X(13:15,:)')
    title('Gyro Biases')
end
%%
outdata = EKF(X_n, P_n,Data,EKF_settings);
% outdata = RTS2(outdata,EKF_settings);
outdata_n = EKF(X, P,Data,EKF_settings);

pos_estimate_n = outdata_n.X(1:3,floor(length(Data.IMU.t)/3*2))';
pos_estimate_conf_n = sqrt(outdata_n.diag_P(1:3,floor(length(Data.IMU.t)/3*2)))';
notsmoothed_position = [pos_estimate_n pos_estimate_conf_n]

pos_estimate = outdata.X(1:3,floor(length(Data.IMU.t)/3*2))';
pos_estimate_conf = sqrt(outdata.diag_P(1:3,floor(length(Data.IMU.t)/3*2)))';
smoothed_position = [pos_estimate pos_estimate_conf]

[~,n] = max(diff(Data.GNSS.t));
drop = [Data.GNSS.pos(:,n)' Data.GNSS.accuracy(1,n) Data.GNSS.accuracy(:,n)']
pick = [Data.GNSS.pos(:,n+1)' Data.GNSS.accuracy(1,n+1) Data.GNSS.accuracy(:,n+1)']

[pos_estimate_n pos_estimate_conf_n; pos_estimate pos_estimate_conf; drop; pick]
%%
figure(10);
sgtitle('Position Result Dataset 3 - Iterative EKF')
% Plot
subplot(3,1,1)
plot(outdata.t,outdata.X(1,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(1,:),'.k')
%plot(Data.NILUS.time(Data_N.NILUS.gps_state),Data.NILUS.position_raw_local(Data_N.NILUS.gps_state,1),'.g')
% hold off 
title('North')
xlabel('time (s)')
ylabel('pos (m)')
subplot(3,1,2)
plot(outdata.t,outdata.X(2,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(2,:),'.k')
%plot(Data.NILUS.time(Data_N.NILUS.gps_state),Data.NILUS.position_raw_local(Data_N.NILUS.gps_state,2),'.g')
% hold off  
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
sgtitle('Position Result Dataset 3 - EKF')
% Plot
subplot(3,1,1)
plot(outdata.t,outdata.X(1,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(1,:),'.k')
%plot(Data.NILUS.time(Data_N.NILUS.gps_state),Data.NILUS.position_raw_local(Data_N.NILUS.gps_state,1),'.g')
% hold off 
title('North')
xlabel('time (s)')
ylabel('pos (m)')
subplot(3,1,2)
plot(outdata.t,outdata.X(2,:))
hold on 
plot(Data.GNSS.t,Data.GNSS.pos(2,:),'.k')
%plot(Data.NILUS.time(Data_N.NILUS.gps_state),Data.NILUS.position_raw_local(Data_N.NILUS.gps_state,2),'.g')
% hold off  
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