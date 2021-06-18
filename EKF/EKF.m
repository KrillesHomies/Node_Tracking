%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code is an Extended Kalman Filter implementation with bias
% estimation. Functions being used in this process can be found in the
% functions folder. 
%
% Input -   initial_state - Is the initial state consiting of 15x1 vector,
%                           with the Position, Velocity, Orientation
%                           (Euler), Accelerometer Bias and Gyroscope Bias
%           initial_Process_covariance - The original uncertainity for each
%                                       state varibale and the covariance,
%                                       the format should be a 15x15 matrix
%           Data - The measurements and inputs with timestamps, see mat
%                  files for examples
%           EKF_settings - All parameters and settings for the filter
%
% Ouput - outdata - Will conatain the predicted state and covariance and
%                   the update state and covariance. 
%
% Edited by Kristoffer Lindve 2021-06-08
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [outdata]= ZUPT_AIDED_EKF(initial_state, initial_Process_covariance,Data,EKF_settings)

% Set intial Matrixes
I = eye(3);
O = zeros(3);
steps = size(Data.IMU.t,1);
% Process noise matrix
Q=blkdiag(EKF_settings.sigma_acc.^2.*I, ...
          EKF_settings.sigma_gyro.^2*I, ...
          EKF_settings.acc_bias_driving_noise.^2*I, ...
          EKF_settings.gyro_bias_driving_noise.^2*I);
      
% Initiate state Values
P = initial_Process_covariance;
X = initial_state;
      
% Allocate Space
[outdata] = allocate_space(X,P,steps);

% Rotation MAtrices
% R_N2B = Rt2b(X(7:9)); %R_N2B = eul2rotm(X(7:9)','XYZ');
% R_B2N = R_N2B';
quat = dcm2q(Rt2b(X(7:9))');

% Input Data
u = [Data.IMU.acc; Data.IMU.gyro]; 

% Setup measurments
GPS.MEAS = Data.GNSS;
GPS.MEAS.COUNT = 0;
ZUPT.MEAS = Data.ZUPT;
ZUPT.MEAS.COUNT = 0;
if(isempty(Data.IMU.mag))
    MAG.MEAS.t = [];
else 
    MAG.MEAS.t = Data.IMU.t;
end
MAG.MEAS.MAG = Data.IMU.mag;
MAG.MEAS.COUNT = 0;
PRESSURE.MEAS = Data.PRESS_SEN;
PRESSURE.MEAS.COUNT = 0;

if(isempty(Data.COMPASS))
    COMPASS.MEAS.t = [];
else 
    COMPASS.MEAS = Data.COMPASS;
end
COMPASS.MEAS.COUNT = 0;

% Motion model
dt = EKF_settings.dt;
g = EKF_settings.g;

% Initialize time at zero
T = Data.IMU.t(1);
latest_pressure = 0;
for k=2:steps
    % Itterate time
    T = T+dt;
    
    % Get input
    u_h = u(:,k-1) + X(10:15);
    outdata.u_h(:, k-1) = u_h;
    %if(~Data.ZUPT.ZUPT(k))
        [X,quat] = nav_eq(X,u_h,quat,g,dt);
    %end
    % Update Covariance
    [F G]=state_matrix(quat,u_h',dt);
    P=F*P*F'+G*Q*G';
    P=(P+P')/2;
    
    % Save Predicted State Variables
    outdata.X_predicted(:, k) = X;
    outdata.P_predicted(:, :, k) = P;
    
    % Get measurement - Check for new values
    H = []; % Linearized model
    Measured = []; % measured Value
    Predicted = []; % predicted Value
    R = []; % Measurment noise 
    % Create Measrment models
    [GPS.MEAS.COUNT GPS.MEAS.NEW_VALUES] = get_next_meas(T,GPS.MEAS.t,GPS.MEAS.COUNT);
    if(GPS.MEAS.NEW_VALUES)
        H_gps = [eye(3) zeros(3,12)];
        Measured = [Measured; Data.GNSS.pos(:,GPS.MEAS.COUNT)];
        %R = [R; GPS.MEAS.pos_acc(:,GPS.MEAS.COUNT).^2];
        acc = Data.GNSS.accuracy(:,GPS.MEAS.COUNT);
        R = blkdiag(R, diag([acc(1) acc(1) acc(2)].^2));
        Predicted = [Predicted; H_gps*X];
        H = [H; H_gps];
    end
    [PRESSURE.MEAS.COUNT PRESSURE.MEAS.NEW_VALUES] = get_next_meas(T,PRESSURE.MEAS.t,PRESSURE.MEAS.COUNT);
    if(PRESSURE.MEAS.NEW_VALUES)
        latest_pressure = PRESSURE.MEAS.press(PRESSURE.MEAS.COUNT);
        if(latest_pressure > EKF_settings.compass_pressure_cutoff)
            
            H_press = [0 0 EKF_settings.pressure*norm(EKF_settings.g)/100 zeros(1,12)];
            Measured = [Measured; PRESSURE.MEAS.press(PRESSURE.MEAS.COUNT)];
            R = blkdiag(R, diag(EKF_settings.sigma_pressure.^2));
            Predicted = [Predicted; H_press*X];
            H = [H; H_press];
        end
    end

    [ZUPT.MEAS.COUNT ZUPT.MEAS.NEW_VALUES] = get_next_meas(T,ZUPT.MEAS.t,ZUPT.MEAS.COUNT);
    if(ZUPT.MEAS.NEW_VALUES)
        if(ZUPT.MEAS.ZUPT(ZUPT.MEAS.COUNT))
            H_zupt = [zeros(3,3) eye(3) zeros(3,9)];
            Measured = [Measured; zeros(3,1)];
            R = blkdiag(R, diag(EKF_settings.sigma_zero_velocity.^2.*ones(3,1)));
            Predicted = [Predicted; H_zupt*X];
            H = [H; H_zupt];
        end
    end
    
  
    [COMPASS.MEAS.COUNT COMPASS.MEAS.NEW_VALUES] = get_next_meas(T,COMPASS.MEAS.t,COMPASS.MEAS.COUNT);
    if(COMPASS.MEAS.NEW_VALUES && latest_pressure > EKF_settings.compass_pressure_cutoff)
        R_N2B = q2dcm(quat)';
        H_mag = [zeros(3,6) eye(3) zeros(3,6)];
        xchi = (H_mag*P*H_mag'+diag(EKF_settings.sigma_mag.^2.*ones(3,1)));
        v = (COMPASS.MEAS.Data(:,COMPASS.MEAS.COUNT)) - X(7:9);
        if(v'*xchi*v < EKF_settings.compass_cutoff) % 90%
            Measured = [Measured; (COMPASS.MEAS.Data(:,COMPASS.MEAS.COUNT))];
            Predicted = [Predicted; X(7:9)];
            R = blkdiag(R, diag(EKF_settings.sigma_mag.*ones(3,1)));
            H = [H; H_mag];
        end
        
    end    

    if(~isempty(R))
        % Calculate the Kalman filter gain.
        K=(P*H')/(H*P*H'+(R));
        % Calculate Update
        z=K*(Measured-Predicted);
        % Correct the navigation states with adjustments
        X=X+z;
        % Correct orientation
        quat=Gamma(quat,z(7:9));
        X(7:9) = dcm2euler(q2dcm(quat)); 

        % Update the Kalman filter state covariance.
        P=(eye(15)-K*H)*P;
        P=(P+P')/2;
    end    
    
    % Save the data to the output data structure
    outdata.X(:,k) = X;
    outdata.P(:,:,k) = P;
    outdata.diag_P(:,k)=diag(P);
    outdata.F(:,:,k)=F;
    outdata.t(k) = T;
    
    
end


function [count,new_values] = get_next_meas(T,t,count)
    new_values = 0;
    if(~isempty(t))
        if(length(t) > count)
            while(round(T,5)>=round(t(count+1),5))
                count = count + 1;
                new_values = 1;
                if(length(t) == count)
                        break;
                end
            end
        end
    end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Error correction of quaternion    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q=Gamma(q,epsilon)

    R = q2dcm(q);
    % Construct skew symetric matrx
    OMEGA=[0 -epsilon(3) epsilon(2); epsilon(3) 0 -epsilon(1); -epsilon(2) epsilon(1) 0];

    % Cortect the DCM matrix
    R=(eye(3)-OMEGA)*R;
    
    q=dcm2q(R);

return

%%
% Allocates space for each of the outputs
%%
function [outdata] = allocate_space(X,P,steps)
    % State - Predicted & Updated 
    outdata.X = zeros(15, steps);
    outdata.X_predicted = zeros(15, steps);
    % Process Covariance - Predicted & Updated
    outdata.P = zeros(15, 15, steps);
    outdata.P_predicted = zeros(15, 15, steps);
    outdata.diag_P = zeros(15, steps);
    % Acceleration Vector once gravity is removed
    outdata.a_N = zeros(6, steps);
    % Measurments 
    outdata.u_h = zeros(6, steps);
    outdata.F = zeros(15, 15, steps);
    % Time vector
    outdata.t = zeros(1, steps);
    % Add initial Values
    outdata.X(:, 1) = X;
    outdata.X_predicted(:, 1) = X;
    outdata.P(:, :, 1) = P;
    outdata.P_predicted(:, :, 1) = P;
    outdata.diag_P(:, 1) = diag(P);
return
