function [outdata_1,outdata_2] = cut_data_imu(Data,time)

imu = Data.IMU;

outdata_1 = Data;
outdata_2 = Data;

[~,i] = min(abs(imu.t - time));

imu_1.t = imu.t(1:i);
imu_1.acc = imu.acc(:,1:i);
imu_1.gyro = imu.gyro(:,1:i);
imu_1.mag = imu.mag(:,1:i);

l = length(imu.t);
imu_2.t = imu.t(i+1:l);
imu_2.acc = imu.acc(:,i+1:l);
imu_2.gyro = imu.gyro(:,i+1:l);
imu_2.mag = imu.mag(:,i+1:l);

outdata_1.IMU = imu_1;
outdata_2.IMU = imu_2;