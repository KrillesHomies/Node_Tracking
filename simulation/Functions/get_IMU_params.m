function [acc_params gyro_params mag_params] = get_IMU_params()
%% Set up IMU 
% Accelerometer
conversion_mg_ms2 = 9.80665;
Acc_Meas_Range = inf; %infinity, no limmit
Acc_Res = 0.02*conversion_mg_ms2; %from data sheet
Acc_const_bias = 0.002*conversion_mg_ms2;
Acc_Noise_den = 20e-9*conversion_mg_ms2;
Acc_bias_instab = 0.04*conversion_mg_ms2;
acc_params = accelparams('MeasurementRange',Acc_Meas_Range, ...
                            'Resolution',Acc_Res, ...
                            'ConstantBias',Acc_const_bias, ...
                            'NoiseDensity',Acc_Noise_den, ...
                            'BiasInstability',Acc_bias_instab);
% Gyro
conversion_deg_rad =pi/180;
Gyro_Meas_Range = inf; %infinity, no limmit
Gyro_Res = 0.003*conversion_deg_rad; %from data sheet
Gyro_bias_instab = 8*conversion_deg_rad/60; % / 3600
Gyro_const_bias = 0.04*conversion_deg_rad;
Gyro_Noise_den = 0.005*conversion_deg_rad;
gyro_params = gyroparams('MeasurementRange',Gyro_Meas_Range ...
                            ,'Resolution',Gyro_Res,...
                            'ConstantBias',Gyro_Noise_den, ...
                            'NoiseDensity',Acc_Noise_den, ...
                            'BiasInstability',Gyro_bias_instab);
% Mag
conversion_guass_utesla = 10e-5/10e-9;
conversion_guass_utesla = 1;
Mag_Meas_Range = inf; %infinity, no limmit
Mag_const_bias = 0.003*conversion_guass_utesla;
Mag_Noise_den = 400*10e-9*conversion_guass_utesla;
mag_params = magparams('MeasurementRange',Mag_Meas_Range, ...
                        'ConstantBias',Mag_const_bias, ...
                        'NoiseDensity',Mag_Noise_den);
                  