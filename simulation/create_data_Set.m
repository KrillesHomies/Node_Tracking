% This script creates a dataset with simulated data to be used in the EKF
% The set is saved in the variable "Dataset"

%% Stettings for simulation
% Dataset
Dataset = []; % Clear Variable
addpath Functions;

% Time step in simulation
dt = 0.01;
% Easy access if you want to change the depth
depth = 40;
%Amount of Datasets create with same initial parameters
size_data_set = 1;
%% Noise settings, changing these will require changes to the system
% Gets the parames for the IMU - To change paramerters edit the function
Sensor = [];
[Sensor.imu.acc_params Sensor.imu.gyro_params Sensor.imu.mag_params] = get_IMU_params();
% load settings, these should also be used for the EKF
EKF_settings = get_sim_noise_settings(); 
%% Set Node attributes
shape = [];
shape.Mass = 50;
shape.Volume = 0.0350;
shape.Height = 0.8;
shape.Radius = 0.4;
shape.Center_of_Bouyancy_sink = [0.05 0.00 0.05]';
shape.Center_of_Bouyancy_rise = [0.05 0.05 0.05]';

shape.Center_of_mass = [-0.05 0.05 0.01]';
shape.Center_of_Pressure = [-0.02 -0.02 -0.02]';
shape.Inertia = blkdiag(3*shape.Mass/20*((2*shape.Radius/3)^2 + 4*shape.Height), ...
                           3*shape.Mass/20*((2*shape.Radius/3)^2 + 4*shape.Height),  ...
                           (2*shape.Radius/3)^2*3/10*shape.Mass );

%% Set Envirmoment parameters
ocean = [];
ocean.depth = depth;
ocean.rotation_cycle_time = 1200;
ocean.amplitude_cycle_time = 500;
ocean.current_rotation_intial = [0 0 pi/4]';
ocean.current_rotation_timevarying = [0 0 pi*0.15]';
ocean.current_rotation_noise = [0 0 0.01]';
ocean.current_amplitude_intial = [0.8 0 0]';
ocean.current_amplitude_timevarying = [0.2 0 0]';
ocean.current_amplitude_noise  = [0.05 0.05 0.05]';

%% Set Constants
Constants = [];
Constants.g = [0 0 9.81];
Constants.p = 1003;
Constants.Timer_bottom = 50;
Constants.Timer_Floating = 5;
% Volume when risning/sinking
Constants.v_sinking = shape.Volume;
Constants.v_rising = 3*shape.Volume;
% Drag coeficent
Constants.C_sinking = 0.4955; 
Constants.C_rising = 0.7955;

%% Set up IMU
Sensor.imu.noNoise = imuSensor('accel-gyro-mag');
Sensor.imu.Noisy = imuSensor('accel-gyro-mag','Accelerometer',Sensor.imu.acc_params,'Gyroscope',Sensor.imu.gyro_params,'Magnetometer',Sensor.imu.mag_params);
Sensor.imu.sampleRate = 1/100;
Sensor.pressuresenor.random_noise = 200;
Sensor.pressuresenor.sampleRate = 1/5;
Sensor.gps.random_noise = 10;
Sensor.gps.sampleRate = 1/4;

%% Set Initial Values
state_i = [];
state_i.A = [0 0 0]';
state_i.Qw = [0 0 0]';
state_i.Qa = [0 0 0]';
state_i.P = EKF_settings.factp(1)*randn(3,1)*0;
state_i.V = EKF_settings.factp(2)*randn(3,1)*0;
state_i.Eul = EKF_settings.factp(3:5)'.*randn(3,1)*0;
state_i.Eul = dcm2euler(Rt2b(state_i.Eul'))'*0;
state_i.Q = dcm2q(Rt2b(state_i.Eul)');


constant_accbias = EKF_settings.factp(6)*randn(3,1);
constant_gyrobias = EKF_settings.factp(7)*randn(3,1);
state_i.accbias = constant_accbias;
state_i.gyrobias  = constant_gyrobias;

%% Create data set
for b=1:size_data_set
    % Initial state
    state = state_i; 
    constant_accbias = EKF_settings.factp(6)*randn(3,1);
    constant_gyrobias = EKF_settings.factp(7)*randn(3,1);
    state_i.accbias = constant_accbias;
    state_i.gyrobias  = constant_gyrobias;
    
    initial_state = [EKF_settings.factp(1)*randn(3,1); ...
                    EKF_settings.factp(2)*randn(3,1); ...
                    dcm2euler(Rt2b(EKF_settings.factp(3:5)'.*randn(3,1)')')'; ... 
                    [0 0 0]';
                    [0 0 0]'];
    % Intialize state
    step_rv = 0;
    % Set time to zero
    T = 0;
    % Initialize values
    timer_value = Constants.Timer_bottom;
    timer2_value = Constants.Timer_Floating;
    Data_N = [];
    Data_NN = [];
    j = 0;
    m = 0;
    n = 0;
    true_state = [];

    for i=1:31998
        % Saves true state
        true_state(:,i) = [state.P; state.V; state.Eul; state.accbias; state.gyrobias];
        % Updates bias
        state.accbias = state.accbias + EKF_settings.acc_bias_driving_noise*randn(3,1);
        state.gyrobias  = state.gyrobias + EKF_settings.gyro_bias_driving_noise*randn(3,1);
        % Selects state
        if(abs(state.P(end)) >=  ocean.depth && step_rv == 0 && timer_value > 0)
            step_rv = 1;
        elseif (timer_value > 0 && step_rv == 1)
            timer_value = timer_value - dt;
        elseif(state.P(end) < 0 && step_rv == 2)
            timer2_value = timer2_value - dt;
            if(timer2_value <= 0)
                break;
            end
        elseif(timer_value <= 0)
            step_rv = 2;
        end

        if(step_rv == 0 || step_rv == 2)
            % Calculate Forces
            [fg, fb, fd, tg, tb ,td] = calculate_forces(state, shape, ocean, Constants,step_rv, T);

            % add some angular velocity  resistans term for stability
            tw = shape.Height*Constants.p*(state.Qw).^2.*sign(state.Qw);
            tw2 = sqrt(abs(state.Qw)).*sign(state.Qw);

            % Calculate Acceleration
            Acc = (fg + fb + fd)./shape.Mass;
            Alpha = diag(inv(shape.Inertia)).*((tg + tb + td) - tw - tw2);
        end

        if(step_rv == 1)
            Alpha = state.Qw/(-dt)*0.5;
            Acc = state.V/(-dt)*0.5;
        end

        if(max(Alpha) > 1000)
            pause(1);
        end
        
        % Round values with resolution
        Acc = round(Acc/Sensor.imu.acc_params.Resolution,1)*Sensor.imu.acc_params.Resolution;
        Alpha = round(Alpha/Sensor.imu.gyro_params.Resolution,1)*Sensor.imu.gyro_params.Resolution;

        % Get IMU data
        [acc_nn, gyro_nn, mag_nn] = Sensor.imu.noNoise(-Acc',-(state.Qw)',q2dcm(state.Q)');
        [acc_n, gyro_n, mag_n] = Sensor.imu.Noisy(-Acc',-(state.Qw)',q2dcm(state.Q)');
%         acc_nn = (q2dcm(state.Q)'*(Acc + EKF_settings.g'))';
%         gyro_nn = (q2dcm(state.Q)'*state.Qw)';
%         mag_nn = (q2dcm(state.Q)'*EKF_settings.MagneticField')';
%         
%         acc_n = acc_nn - state.accbias' + EKF_settings.sigma_acc*randn(1,3);
%         gyro_n = gyro_nn - state.gyrobias' + EKF_settings.sigma_gyro*randn(1,3);
%         mag_n = mag_nn  + EKF_settings.sigma_mag*randn(1,3);
        acc_n = acc_nn - state.accbias';
        gyro_n = gyro_nn - state.gyrobias';
        
    
        x = [state.P; state.V; state.Eul; state.accbias; state.gyrobias];
        [x,state.Q] = nav_eq(x,[acc_nn gyro_nn]',state.Q,EKF_settings.g,dt);
        
%         A_navigation_frame = R_B2N*acc_nn' - Constants.g';
%         
%         % Navigation Algorithm 
%         A = [1 0 0 dt 0 0; ...
%              0 1 0 0 dt 0; ...
%              0 0 1 0 0 dt; ...
%              0 0 0 1 0 0; ...
%              0 0 0 0 1 0; ...
%              0 0 0 0 0 1];
%         
%         B = [0.5*dt^2 0        0; ...
%              0        0.5*dt^2 0; ...
%              0        0        0.5*dt^2; ...
%              dt       0        0; ...
%              0        dt       0; ...
%              0        0        dt];
%          
%         x = A*[state.P; state.V] + B*Acc;
% 
%         P=gyro_nn(1)*dt;
%         Q=gyro_nn(2)*dt;
%         R=gyro_nn(3)*dt;
% 
%         OMEGA=zeros(4);
%         OMEGA(1,1:4)=0.5*[0 R -Q P];
%         OMEGA(2,1:4)=0.5*[-R 0 P Q];
%         OMEGA(3,1:4)=0.5*[Q -P 0 R];
%         OMEGA(4,1:4)=0.5*[-P -Q -R 0];
% 
%         v=norm(gyro_nn)*dt;
%         state.Q = dcm2q(R_B2N);
%         if v~=0
%             state.Q=(cos(v/2)*eye(4)+2/v*sin(v/2)*OMEGA )*state.Q;
%         end
%         
%         R_B2N = q2dcm(state.Q);
%         R_N2B = R_B2N';
        
        
        % Pressure Sesnors without noise
        p = abs(state.P(end))*Constants.p*norm(Constants.g)/100;
        if(p < 10 || state.P(end) < 0)
            p = 1003/100;
        end

        % Update Angle
        state.P = x(1:3);
        state.V = x(4:6);
        state.Eul = dcm2euler(q2dcm(state.Q))'; %state.Eul = rotm2eul(R_N2B,'XYZ')';
        state.Qw = state.Qw + Alpha*dt;

        T = round(T + dt,5);

        % Save Data
        if (mod(T, Sensor.gps.sampleRate) == 0)
            j = j+1;
            Data_NN.GNSS.pos(:,j) = state.P';
            Data_NN.GNSS.t(j,:) = T;
            Data_NN.GNSS.accuracy(:,j) = [EKF_settings.sigma_gps(1) EKF_settings.sigma_gps(2)];
            Data_N.GNSS.pos(:,j) = state.P' + EKF_settings.sigma_gps.*randn(1,3);
            Data_N.GNSS.t(j,:) = T;
            Data_N.GNSS.accuracy(:,j) = [EKF_settings.sigma_gps(1) EKF_settings.sigma_gps(2)];
        end
        if (mod(T, Sensor.imu.sampleRate) == 0)
            m  = m +1;
            Data_NN.IMU.t(m,:) = T;
            Data_NN.IMU.acc(:,m) = acc_nn;
            Data_NN.IMU.gyro(:,m) = gyro_nn;
            Data_NN.IMU.mag(:,m) = mag_nn;
            Data_N.IMU.t(m,:) = T;
            Data_N.IMU.acc(:,m) = acc_n;
            Data_N.IMU.gyro(:,m) = gyro_n;
            Data_N.IMU.mag(:,m) = mag_n;
        end
        if (mod(T, Sensor.pressuresenor.sampleRate) == 0)
            n = n+1;
            Data_NN.PRESS_SEN.t(n,:) = T;
            Data_NN.PRESS_SEN.press(:,n) = p;
            Data_N.PRESS_SEN.t(n,:) = T;
            Data_N.PRESS_SEN.press(:,n) = p + EKF_settings.sigma_pressure*randn(1,1);
        end
    end
    
    % Add run to dataset array
    Dataset{b}.Data_NN = Data_NN;
    Dataset{b}.Data_N = Data_N;
    Dataset{b}.X_r = true_state;
    Dataset{b}.initial_state = initial_state;
    b
end

disp('Data set is saved in the variable "Dataset", save it where you want.')
