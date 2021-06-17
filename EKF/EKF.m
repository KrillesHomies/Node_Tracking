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
            %R = [R; EKF_settings.sigma_zero_velocity.^2.*ones(3,1)];
            R = blkdiag(R, diag(EKF_settings.sigma_zero_velocity.^2.*ones(3,1)));
            Predicted = [Predicted; H_zupt*X];
            H = [H; H_zupt];
        end
    end
    
    
    %[MAG.MEAS.COUNT MAG.MEAS.NEW_VALUES] = get_next_meas(T,MAG.MEAS.t,MAG.MEAS.COUNT);
    %window_size = 11;
%     % MAG.MEAS.NEW_VALUES && GRAVITY.MEAS.NEW_VALUES &&  && latest_pressure > 100
%     if(k > window_size && steps > k+window_size && ~isempty(Data.IMU.mag) && latest_pressure > 100)
% 
%         accel_data = Data.IMU.acc(:,k-floor(window_size/2):k+ceil(window_size/2)-1,:);
%         mag_data = Data.IMU.mag(:,k-floor(window_size/2):k+ceil(window_size/2)-1,:);
%         V1 = [EKF_settings.g'/norm(EKF_settings.g).*ones(3,window_size)...
%               EKF_settings.MagneticField'/norm(EKF_settings.MagneticField').*ones(3,window_size)];
%         V2 = [accel_data ...
%               mag_data];
%         V2 = V2./sqrt(sum(V2.^2));
%         ROT = V2*V1';
%         [U S V] = svd(ROT);
%         corrections = flip(rotm2eul(V*U'));  
%         H_mag_roll_pitch = [zeros(2,6) eye(2) zeros(2,7)];
%         H_mag_yaw = [zeros(1,8) eye(1) zeros(1,6)];
%         v = (corrections-X(7:9)');
%         xchi_rp = (H_mag_roll_pitch*P*H_mag_roll_pitch'+diag(EKF_settings.sigma_mag.^2.*ones(2,1)));
%         xchi_y = (H_mag_yaw*P*H_mag_yaw'+diag(EKF_settings.sigma_mag.^2.*ones(1,1)));
%         g_mag_test = abs(1-mean(sqrt(sum(accel_data.^2)))/norm(EKF_settings.g'));
%         if(g_mag_test < 0.05)
%             if(v(1:2)*xchi_rp*v(1:2)' < 10e-5) %10e-5
%                 Measured = [Measured; X(7:8)];
%                 Predicted = [Predicted; corrections(1:2)'];
%                 R = blkdiag(R, diag(EKF_settings.sigma_mag.^2.*ones(2,1)));
%                 H = [H; H_mag_roll_pitch];
%             end
%             if( v(3)*xchi_y*v(3)' < 10e-5) %10e-5 0.58
%                 Measured = [Measured; X(9)];
%                 Predicted = [Predicted; corrections(3)];
%                 R = blkdiag(R, diag(EKF_settings.sigma_mag.^2.*ones(1,1)));
%                 H = [H; H_mag_yaw];
%             end
%         end
%     end 
  
    [COMPASS.MEAS.COUNT COMPASS.MEAS.NEW_VALUES] = get_next_meas(T,COMPASS.MEAS.t,COMPASS.MEAS.COUNT);
    if(COMPASS.MEAS.NEW_VALUES && latest_pressure > EKF_settings.compass_pressure_cutoff)
        R_N2B = q2dcm(quat)';
        H_mag = [zeros(3,6) eye(3) zeros(3,6)];
        %H_mag = [zeros(3,6) skew(R_N2B*EKF_settings.MagneticField') zeros(3,6)];
        xchi = (H_mag*P*H_mag'+diag(EKF_settings.sigma_mag.^2.*ones(3,1)));
        v = (COMPASS.MEAS.Data(:,COMPASS.MEAS.COUNT)) - X(7:9);
%         v1 = dcm2q(Rt2b(dcm2euler(quat2rotm(COMPASS.MEAS.Data(COMPASS.MEAS.COUNT,:))))');
%         v2 = dcm2q(Rt2b(X(7:9))');
%         v = dcm2euler(q2dcm(quatmultiply(conj(v1'),v2'))')';
        if(v'*xchi*v < EKF_settings.compass_cutoff) % 90%
            Measured = [Measured; (COMPASS.MEAS.Data(:,COMPASS.MEAS.COUNT))];
            Predicted = [Predicted; X(7:9)];
            R = blkdiag(R, diag(EKF_settings.sigma_mag.*ones(3,1)));
            H = [H; H_mag];
        end
%         H_mag_roll_pitch = [zeros(2,6) eye(2) zeros(2,7)];
%         H_mag_yaw = [zeros(1,8) eye(1) zeros(1,6)];
%         xchi_rp = (H_mag_roll_pitch*P*H_mag_roll_pitch'+diag(EKF_settings.sigma_mag.^2.*ones(2,1)));
%         xchi_y = (H_mag_yaw*P*H_mag_yaw'+diag(EKF_settings.sigma_mag.^2.*ones(1,1)));        
%         if(v(1:2)'*xchi_rp*v(1:2) < 0.001) %10e-5
%             Measured = [Measured; (COMPASS.MEAS.Data(1:2,COMPASS.MEAS.COUNT))];
%             Predicted = [Predicted; X(7:8)];
%             R = blkdiag(R, diag(EKF_settings.sigma_mag.^2.*ones(2,1)));
%             H = [H; H_mag_roll_pitch];
%         end
%         if( v(3)'*xchi_y*v(3) < 0.001) %10e-5 0.58
%             Measured = [Measured; (COMPASS.MEAS.Data(3,COMPASS.MEAS.COUNT))];
%             Predicted = [Predicted; X(9)];
%             R = blkdiag(R, diag(EKF_settings.sigma_mag.^2.*ones(1,1)));
%             H = [H; H_mag_yaw];
%         end
        
    end    
%     [MAG.MEAS.COUNT MAG.MEAS.NEW_VALUES] = get_next_meas(T,MAG.MEAS.t,MAG.MEAS.COUNT);
%     if(MAG.MEAS.NEW_VALUES)
%         R_N2B = q2dcm(quat)';
%         [~, St_1] = Lineariztion_MAG(dcm2euler(q2dcm(quat)), EKF_settings.MagneticField');
%         H_mag = [zeros(3,6) St_1 zeros(3,6)];
%         %H_mag = [zeros(3,6) skew(R_N2B*EKF_settings.MagneticField') zeros(3,6)];
%         xchi = (H_mag*P*H_mag'+diag(EKF_settings.sigma_mag.^2.*ones(3,1)));
%         v = (MAG.MEAS.MAG(:,MAG.MEAS.COUNT)-R_N2B*EKF_settings.MagneticField');
%         %if(v'*xchi*v < 0.58) % 90%
%             Measured = [Measured; MAG.MEAS.MAG(:,MAG.MEAS.COUNT) ];
%             Predicted = [Predicted; R_N2B*EKF_settings.MagneticField'];
%             R = blkdiag(R, diag(EKF_settings.sigma_mag.*ones(3,1)));
%             H = [H; H_mag];
%         %end
%     end
%     [GRAVITY.MEAS.COUNT GRAVITY.MEAS.NEW_VALUES] = get_next_meas(T,GRAVITY.MEAS.t,GRAVITY.MEAS.COUNT);
%     if(GRAVITY.MEAS.NEW_VALUES)
%         if(GRAVITY.MEAS.MAG(GRAVITY.MEAS.COUNT))
%             R_N2B = q2dcm(quat)';
%             [~, St_2] = Lineariztion_MAG(dcm2euler(q2dcm(quat)), EKF_settings.g');
%             H_mag = [zeros(3,6) St_2 eye(3)*0 zeros(3,3)];
%             xchi = (H_mag*P*H_mag'+diag(EKF_settings.sigma_acc.^2.*ones(3,1)));
%             v = (u_h(1:3)-(R_N2B*(-EKF_settings.g')));
%             if(v'*xchi*v < 0.58) % 90%
%                 Measured = [Measured; u_h(1:3)];
%                 Predicted = [Predicted; (R_N2B*(-EKF_settings.g'))];
%                 R = blkdiag(R, diag(EKF_settings.sigma_acc.*ones(3,1)));
%                  H = [H; H_mag];
%             end
%         end
%     end

    if(~isempty(R))
        % Calculate the Kalman filter gain.
        K=(P*H')/(H*P*H'+(R));
        % Update the perturbation state estimate.
        z=K*(Measured-Predicted);
        % Correct the navigation states using current perturbation estimates.
        X=X+z;
        

        quat=Gamma(quat,z(7:9));
        
        X(7:9) = dcm2euler(q2dcm(quat)); %X(7:9) = rotm2eul(R_N2B,'XYZ');

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

% function Hn=nummericalJacobian(R_N2B, Vector)
% ori = rotm2eul(R_N2B,'XYZ');
% h = 0.000001;
% Hn=zeros(3);
% for ii=1:3
%     euler=ori;
%     euler(ii)=euler(ii)+h;
%     Hn(:,ii)=(eul2rotm(euler,'XYZ')'-R_N2B')*(Vector')./h;
% end
% return