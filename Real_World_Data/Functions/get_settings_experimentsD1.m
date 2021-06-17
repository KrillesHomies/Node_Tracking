
function settings=get_settings()


% Process noise covariance (Q)
settings.sigma_acc = 0.01; % [m/s^2]
settings.sigma_gyro =0.05*pi/180; % [rad/s]
settings.sigma_acc = 0.01; % [m/s^2]
settings.sigma_gyro =0.01*pi/180; % [rad/s]

% Process noise for modeling the drift in accelerometer biases (x,y,z 
% platform coordinate axis) [m/s^2].
settings.acc_bias_driving_noise=0.0000000001; 
% Process noise for modeling the drift in gyroscope biases (x,y,z platform
% coordinate axis) [rad/s].
settings.gyro_bias_driving_noise=0.0000000001*pi/180; 

% process noise
settings.factp = [1.0 ... % Initial process noise for Position [m]
                 1.0 ... % Initial process noise for Velocity [m/s]
                  (pi/180*[5 5 15]) ... % Initial Aguglar (roll,pitch,yaw) uncertainty [rad], Yaw is higher since it is more uncertain
                  0.1 ... % Accelerometer biases [m/s^2]
                  (0.1*pi/180)]; % Gyro biases [rad/s]  


% Measurment noise               
settings.sigma_gps = sqrt([1.7 1.7 2]); % Measurment noise gps
settings.sigma_pressure = 10.0; % Measurment noise pressure sensor
settings.sigma_mag = 0.5; % Measurment noise Mag
settings.sigma_zero_velocity = 0.00001;

% Zero Velocity Update
settings.window_size = 100;
settings.gamma = 8500;

% Global Enviromental Factors
settings.MagneticField = [27.5550   -2.4169  -16.0849];
settings.g = [0 0 9.81];
settings.pressure = 1005;

% Sampling Rate
settings.dt = 0.01;

end




