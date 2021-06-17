% use gravity vector?
close all
addpath EKF
addpath simulation_data
load('Data_set40.mat')
EKF_settings = get_sim_noise_settings();
% Get Settings 
I = eye(3); 
X_initial = [0 0 0 ... % Initial Position
                 0 0 0 ... % Initial Velocity
                 0 0 0 ... % Initial Attitude (euler)
                 0 0 0 ... % Initial Acc Bias
                 0 0 0]';  % Initial Gyro Bias 
P_initial = blkdiag(EKF_settings.factp(1)^2*I, ...
            EKF_settings.factp(2)^2*I, ...
            diag(EKF_settings.factp(3:5)).^2, ...
            EKF_settings.factp(6)^2*I, ...
            EKF_settings.factp(7)^2*I);  

X = X_initial; 
P = P_initial;
X_n = X;
P_n = P;

Data_vector_EKF = [];
Data_vector_smoothed = [];
pos_diff = [];
pos_diff_s = [];
t_half = floor(length(Dataset{1}.Data_N.IMU.t)/2);
t_half_gnss = floor(length(Dataset{1}.Data_NN.GNSS.t(:))/2);

for j=1:size(Dataset,2)

    Data = Dataset{j}.Data_N;
    X = Dataset{j}.initial_state;
    %% Remove GPS points which are under the water
    i = (Dataset{j}.Data_NN.GNSS.pos(3,:) < 0);
    Data.GNSS.pos = Data.GNSS.pos(:,i);
    Data.GNSS.t = Data.GNSS.t(i);
    %% ZUPT
    EKF_settings = get_sim_noise_settings();
    EKF_settings.gamma = 5;
    EKF_settings.window_size = 10;
    [T_zupt, zupt] = ZUPT_CAL([Data.IMU.acc; Data.IMU.gyro],EKF_settings,'GLRT');

    % Add to Data
    Data.ZUPT.t = Data.IMU.t;
    Data.ZUPT.ZUPT = zupt;
    %% Compass
    Data.COMPASS.t = [];
    Data.COMPASS.data = [];
    Data.COMPASS = get_compass(Data.IMU.acc,Data.IMU.mag,Data.IMU.t,11,EKF_settings);
    
    %% Settings
    EKF_settings.sigma_acc = 0.50 ; % 60 0.20  0.50
    EKF_settings.sigma_gyro = 0.20*pi/180 ;  % 34 0.10*pi/180 0.20*pi/180
    EKF_settings.sigma_mag = 0.25;  
    EKF_settings.compass_cutoff = 0.0001;
    EKF_settings.compass_pressure_cutoff = 50;
    
    %% Run filter
    Data_vector_EKF{j} = EKF(X, P,Data,EKF_settings);
    Data_vector_smoothed{j} = RTS2(Data_vector_EKF{j},EKF_settings);
    
    %% Set True path
    Data_true = Dataset{j}.X_r;
    data_real = Dataset{j}.Data_NN;
    
    %% Plot
    figure(10);
    sgtitle('EKF')
    % Plot
    subplot(3,1,1)
    plot(data_real.GNSS.t(:),data_real.GNSS.pos(1,:),'.k')
    title('North')
    hold on
    plot(Data_vector_EKF{j}.t,Data_vector_EKF{j}.X(1,:))
    ylabel('m')
    xlabel('t(s)')
    subplot(3,1,2)
    plot(data_real.GNSS.t(:),data_real.GNSS.pos(2,:),'.k')
    title('East')
    hold on
    plot(Data_vector_EKF{j}.t,Data_vector_EKF{j}.X(2,:))
    ylabel('m')
    xlabel('t(s)')
    subplot(3,1,3)
    plot(data_real.GNSS.t(:),data_real.GNSS.pos(3,:),'.k')
    title('Down')
    hold on
    plot(Data_vector_EKF{j}.t,Data_vector_EKF{j}.X(3,:))
    ylabel('m')
    xlabel('t(s)')
    

    figure(12);
    scatter(Dataset{j}.X_r(1,t_half),Dataset{j}.X_r(2,t_half),'*k');
    hold on 
    scatter(Data_vector_EKF{j}.X(1,t_half),Data_vector_EKF{j}.X(2,t_half),'MarkerFaceColor',colors(:,j));
    e = ellipse(Data_vector_EKF{j}.X(1:2,t_half)',sqrt(Data_vector_EKF{j}.P(1,1,t_half)),sqrt(Data_vector_EKF{j}.P(2,2,t_half))); 
    plot((e(:,1)),(e(:,2)),'-.','color',colors(:,j));
    xlabel('North')
    ylabel('East')
    title('EKF model with 1 uncertainty bound')
    
    pos_diff(j) = norm(Dataset{j}.X_r(1:3,t_half) - Data_vector_EKF{j}.X(1:3,t_half));
    pos_diff_s(j) = norm(Dataset{j}.X_r(1:3,t_half) - Data_vector_smoothed{j}.X(1:3,t_half));

end

T = [];
T.avg = mean(pos_diff);
T.std = std(pos_diff);
T.max = max(pos_diff);
T.min = min(pos_diff);
T.Distxy = norm(Dataset{j}.X_r(1:2,t_half));
T.Distxyz = norm(Dataset{j}.X_r(1:3,t_half));

T_s = [];
T_s.avg = mean(pos_diff_s);
T_s.std = std(pos_diff_s);
T_s.max = max(pos_diff_s);
T_s.min = min(pos_diff_s);


% T_s
%[T_s.avg T_s.std T_s.max T_s.min T.Distxy T.Distxyz]
[T.avg T.std T.max T.min T.Distxy T.Distxyz]


