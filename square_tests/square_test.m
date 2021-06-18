%% Runs code on test data of a square
% The Final and initial position should be about the same.
% Two squares are drawn with the IMU
% Run code to see results

close all
addpath(genpath('EKF'))
addpath(genpath('simulation'))
addpath(genpath('square_tests'))
load('square1.mat')
%load('square2.mat')
%% Get Settings 
I = eye(3);
EKF_settings = get_settings_square();
initial_state = [0 0 0 ... % Initial Position
                 0 0 0 ... % Initial Velocity
                 dcm2euler(quat2rotm(quat(:,1)')) ... % Initial Attitude (euler)
                 0 0 0 ... % Initial Acc Bias
                 0 0 0]';  % Initial Gyro Bias
P_initial = blkdiag(EKF_settings.factp(1)^2*I, ...
            EKF_settings.factp(2)^2*I, ...
            diag(EKF_settings.factp(3:5)).^2, ...
            EKF_settings.factp(6)^2*I, ...
            EKF_settings.factp(7)^2*I);  
X = initial_state; 
P = P_initial; 
clear Data
Data.IMU.t = (1:length(u(1:3,500:end)))'*0.01;
Data.IMU.acc = u(1:3,1:end);
Data.IMU.gyro = u(4:6,1:end);
Data.IMU.mag = [];
%% Initial ORientation
% Under the assumption that the system is stationary during the first 20
% samples, the initial roll and pitch is calculate from the 20 first
% accelerometer readings.
f_u=mean(u(1,1:20));
f_v=mean(u(2,1:20));
f_w=mean(u(3,1:20));
roll=atan2(-f_v,-f_w);
pitch=atan2(f_u,sqrt(f_v^2+f_w^2));
% Set the attitude vector
attitude=[roll pitch 0]';
% Calculate quaternion corresponing to the initial attitude
Rb2t=Rt2b(attitude)';
quat=dcm2q(Rb2t);
initial_state(7:9) = attitude;
%% Intiate empty data for filter, these will be ignored in the filter
Data.GNSS = [];
Data.GNSS.t = [];
Data.PRESS_SEN = [];
Data.PRESS_SEN.t = [];
Data.COMPASS.t = [];
Data.COMPASS.MAG = [];

%% ZUPT
EKF_settings.window_size = 10;
EKF_settings.gamma = 25;
[T, zupt] = ZUPT_CAL([Data.IMU.acc; Data.IMU.gyro],EKF_settings,'GLRT');
Data.ZUPT.t = Data.IMU.t;
Data.ZUPT.ZUPT = zupt;

EKF_settings.g = -[0 0 9.8173]; 
%%

% diag(P)'
outdata_ekf = EKF(initial_state, P_initial,Data,EKF_settings);
outdata_ekf = RTS2(outdata_ekf,EKF_settings);

outdata_itt_ekf = outdata_ekf;
for i=1:25
    X = [initial_state(1:9);  outdata_itt_ekf.X(10:end,1)]; 
    P = blkdiag((P_initial(1:9,1:9)), ((outdata_itt_ekf.P(10:end,10:end,end))));
    outdata_itt_ekf = EKF(X, P,Data,EKF_settings);
    outdata_itt_ekf = RTS2(outdata_itt_ekf,EKF_settings); 
   
end

figure(10);
sgtitle('Smoothed Trajectory - EKF')
% Plot
subplot(3,1,1)
plot(outdata_ekf.t,outdata_ekf.X(1,:))
subplot(3,1,2)
plot(outdata_ekf.t,outdata_ekf.X(2,:))
subplot(3,1,3)
plot(outdata_ekf.t,outdata_ekf.X(3,:))

figure(11);
sgtitle('Smoothed Trajectory - Iterative EKF')
% Plot
subplot(3,1,1)
plot(outdata_itt_ekf.t,outdata_itt_ekf.X(1,:))
subplot(3,1,2)
plot(outdata_itt_ekf.t,outdata_itt_ekf.X(2,:))
subplot(3,1,3)
plot(outdata_itt_ekf.t,outdata_itt_ekf.X(3,:))

figure(12);
sgtitle('Smoothed Trajectory - EKF')
plot(outdata_ekf.X(1,:),outdata_ekf.X(2,:))

figure(13);
sgtitle('Smoothed Trajectory - Iterative EKF')
plot(outdata_itt_ekf.X(1,:),outdata_itt_ekf.X(2,:))

display('Desired final position: 0,0');
display(['EKF Final Position: ' num2str(outdata_ekf.X(1,end)) ',' num2str(outdata_ekf.X(2,end))]);
display('Desired final position: 0,0');
display(['Iterative EKF Final Position: ' num2str(outdata_itt_ekf.X(1,end)) ',' num2str(outdata_itt_ekf.X(2,end))]);