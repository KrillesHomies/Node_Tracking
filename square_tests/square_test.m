
close all
% addpath Smoother
addpath EKF
addpath simulation_data
load('square1.mat')
% Get Settings 
I = eye(3);
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

Data.GNSS = [];
Data.GNSS.t = [];
Data.PRESS_SEN = [];
Data.PRESS_SEN.t = [];

Data.COMPASS.t = [];
Data.COMPASS.MAG = [];
EKF_settings = get_settings_square();
EKF_settings.window_size = 10;
EKF_settings.gamma = 25;
[T, zupt] = ZUPT_CAL([Data.IMU.acc; Data.IMU.gyro],EKF_settings,'GLRT');
Data.ZUPT.t = Data.IMU.t;
Data.ZUPT.ZUPT = zupt;

% initial_state(13:15) = mean(Data.IMU.gyro(:,1:100),2);
EKF_settings.g = -[0 0 9.8173]; 
%%

% diag(P)'
outdata = EKF(initial_state, P,Data,EKF_settings);
outdata = RTS2(outdata,EKF_settings);
% initial_state = outdata.X(:,1);
% P = outdata.P(:,:,1);

initial_state(10:15) = outdata.X(10:15,1)
d = diag(outdata.P(:,:,1));
P = diag([diag(P(1:9,1:9)); d(10:15)]);
% initial_state(7:15) = outdata.X(7:15,1)
% d = diag(outdata.P(:,:,1));
% P = diag([diag(P(1:6,1:6)); d(7:15)]);

figure(10);
sgtitle('Position')
% Plot
subplot(3,1,1)
plot(outdata.t,outdata.X(1,:))
subplot(3,1,2)
plot(outdata.t,outdata.X(2,:))
subplot(3,1,3)
plot(outdata.t,outdata.X(3,:))

figure(12);
sgtitle('Orientation')
% Plot
subplot(3,1,1)
plot(outdata.t,outdata.X(7,:))
subplot(3,1,2)
plot(outdata.t,outdata.X(8,:))
subplot(3,1,3)
plot(outdata.t,outdata.X(9,:))


figure(11);
plot(outdata.X(1,:),outdata.X(2,:))
%     
