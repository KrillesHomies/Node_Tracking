%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses singular value decomposition to calculate the orientation using the
% magnetometer data, accelerometer data compared against the expected
% gravity and the expected magnetic vector. The average value over a window
% is used for more stability in the estimate.
%
% Input -   acc - The accelerometer data
%           mag - Magnetometer data
%           t - timestamp for accelerometer and Magnetometer data
%           window_size - size of the window per measurment
%           EKF_settings - All parameters and settings for the filter
%
% Ouput - comp - Orientation estimate in Euler angles with timestamp 
% Edited by Kristoffer Lindve 2021-06-08
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [comp] = get_compass(acc,mag,t,window_size,EKF_settings)

comp.t = [];
comp.Data = [];
i = 1;
for k=1:length(t)
    if(k > window_size && length(t) > k+window_size)
        accel_data = acc(:,k-floor(window_size/2):k+ceil(window_size/2)-1,:);
        mag_data = mag(:,k-floor(window_size/2):k+ceil(window_size/2)-1,:);
        V1 = [EKF_settings.g'/norm(EKF_settings.g).*ones(3,window_size)...
              EKF_settings.MagneticField'/norm(EKF_settings.MagneticField').*ones(3,window_size)];
        V1 = [V1 cross(V1(:,window_size+1:window_size*2),V1(:,1:window_size))];
        V2 = [accel_data ...
              mag_data];
        V2 = [V2 cross(V2(:,window_size+1:window_size*2),V2(:,1:window_size))];  
        V2 = V2./sqrt(sum(V2.^2));
        ROT = V2*V1';
        [U S V] = svd(ROT);
        comp.Data(:,i) = flip(rotm2eul(V*U'))';  
        comp.t(:,i) = t(k);
        i = i + 1;
    end
end