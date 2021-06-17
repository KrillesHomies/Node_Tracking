% use gravity vector?
close all
addpath EKF
addpath simulation_datag
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

% Create Randomized color scheme
colors = linspace(0.1,1.0,size(Dataset,2)).*[1; 1; 1];
colors(1:prod(size(colors))) = colors(randperm(prod(size(colors))));

ellipse = @(point, n_unc, e_unc) [n_unc*cos(-0.1:0.1:(2*pi)); e_unc*sin(-0.1:0.1:(2*pi))]'+point;


for j=1:size(Dataset,2)
%for j=1:5
    % Preprocess
    %j
    Data = Dataset{j}.Data_N;
    X = Dataset{j}.initial_state;
    % Remove GPS points
    i = (Dataset{j}.Data_NN.GNSS.pos(3,:) < 0);
    Data.GNSS.pos = Data.GNSS.pos(:,i);
    Data.GNSS.t = Data.GNSS.t(i);
    % Get ZUPT
    EKF_settings = get_sim_noise_settings();
    EKF_settings.gamma =5;
    EKF_settings.window_size = 10;
    [T_zupt, zupt] = ZUPT_CAL([Data.IMU.acc; Data.IMU.gyro],EKF_settings,'GLRT');

    % Add to Data
    Data.ZUPT.t = Data.IMU.t;
    Data.ZUPT.ZUPT = zupt;

    Data.COMPASS.t = [];
    Data.COMPASS.data = [];
    Data.COMPASS = get_compass(Data.IMU.acc,Data.IMU.mag,Data.IMU.t,11,EKF_settings);

%     EKF_settings.sigma_acc = 0.30; % 60
%     EKF_settings.sigma_gyro = 0.30*pi/180;  % 34
%     EKF_settings.sigma_mag = 0.25;

    EKF_settings.sigma_acc = 0.50 ; % 60 0.20 
    EKF_settings.sigma_gyro = 0.20*pi/180 ;  % 34 0.10*pi/180
    EKF_settings.sigma_mag = 0.25;  
    EKF_settings.compass_cutoff = 0.0001;
    EKF_settings.compass_pressure_cutoff = 50;
    
    X_n = X;
    P_n = P;
    % Smoothing techique
    for o=1:10

        DataS = Data;
        DataS.IMU.t = DataS.IMU.t(1:floor(length(DataS.IMU.t)/2));
        outdata= EKF(X_n, P_n,DataS,EKF_settings);
        soutdata = RTS2(outdata,EKF_settings);

        X_n = (soutdata.X(:,1)+X_n)/(2); 
        P_n = ((((((soutdata.P(:,:,1))))) + (P_n))/(2));
        X_n = (soutdata.X(:,1));  % better sol
        P_n = soutdata.P(:,:,1);
%         X_n = [X(1:6);  X_n(7:end,1)]; 
%         P_n = blkdiag((P(1:6,1:6)), ((P_n(7:end,7:end))));
        X_n = [X(1:9);  X_n(10:end,1)]; 
        P_n = blkdiag((P(1:9,1:9)), ((P_n(10:end,10:end))));
    
    end
    % Run filter
    Data_vector_EKF{j} = EKF(X_n, P_n,Data,EKF_settings);
    Data_vector_smoothed{j} = RTS2(Data_vector_EKF{j},EKF_settings);
    % Set True path
    Data_true = Dataset{j}.X_r;
    data_real = Dataset{j}.Data_NN;
    
    figure(10);
    sgtitle('Iterative EKF')
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
    title('Iterative model with 1 uncertainty bound')
    
    
    pos_diff(j) = norm(Dataset{j}.X_r(1:3,t_half) - Data_vector_EKF{j}.X(1:3,t_half));
    pos_diff_s(j) = norm(Dataset{j}.X_r(1:3,t_half) - Data_vector_smoothed{j}.X(1:3,t_half));

end

T = [];
T.avg = mean(pos_diff);
T.std = std(pos_diff);
T.max = max(pos_diff);
T.min = min(pos_diff);
T.Distxy = norm(Dataset{10}.X_r(1:2,t_half));
T.Distxyz = norm(Dataset{10}.X_r(1:3,t_half));

[T.avg T.std T.max T.min T.Distxy T.Distxyz]



