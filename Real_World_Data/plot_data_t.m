function plot_data(data)

direction={'north','east','down'};
euler_axis={'roll','pitch','yaw'};

figure(1)
clf
h(1)=plot3(data.X(1,:),data.X(2,:),data.X(3,:),'r');
axis equal
grid minor
xlabel('North [m]')
ylabel('East [m]')
zlabel('Downh [m]')
title('Position')
legend(h,{'estimated'})


figure(2)
clf
h=zeros(1,2);
for ii=1:3
subplot(3,1,ii)
hold on;
h(1)=plot(data.t,data.X(ii,:),'r');
h(2)=plot(data.t,data.X(ii,:)+3*sqrt(data.diag_P(ii,:)),'r:');
h(2)=plot(data.t,data.X(ii,:)-3*sqrt(data.diag_P(ii,:)),'r:');
legend(h,{'estimated','3 \sigma'})
ylabel('m')
xlabel('t [s]')
title(['Position ' direction(ii)])
grid minor;
end




figure(3)
clf
h=zeros(1,2);
for ii=1:3
subplot(3,1,ii)
hold on;
h(1)=plot(data.t,data.X(3+ii,:),'r');
h(2)=plot(data.t,data.X(3+ii,:)+3*sqrt(data.diag_P(3+ii,:)),'r:');
h(2)=plot(data.t,data.X(3+ii,:)-3*sqrt(data.diag_P(3+ii,:)),'r:');
legend(h,{'estimated','3 \sigma'})
ylabel('m/s')
xlabel('t [s]')
title(['Speed ' direction(ii)])
grid minor;
end
xlabel('t [s]')

figure(4)
clf
% h=zeros(1,3);
% euler_ang=zeros(3,length(data.t));
% euler_ang_h=zeros(3,length(data.t));
% for jj=1:length(data.t)
%    euler_ang(:,jj)=180/pi*((data.x(7:9,jj))); 
%    euler_ang_h(:,jj)=180/pi*data.X(7:9,jj); 
% end

euler_ang_h =180/pi*data.X(7:9,:); 


for ii=1:3
subplot(3,1,ii)
hold on;
h(1)=plot(data.t, euler_ang_h(ii,:),'r');
h(2)=plot(data.t,euler_ang_h(ii,:)+3*180/pi*sqrt(data.diag_P(6+ii,:)),'r:');
h(2)=plot(data.t,euler_ang_h(ii,:)-3*180/pi*sqrt(data.diag_P(6+ii,:)),'r:');
legend(h,{'estimated','3 \sigma'})
ylabel('deg')
xlabel('t [s]')
title(['Attitude: ' euler_axis(ii)])
grid minor;
end



figure(5)
clf
h=zeros(1,2);
for ii=1:3
subplot(3,1,ii)
hold on;
h(1)=plot(data.t, data.X(ii+9,:),'r');
h(1)=plot(data.t, data.X(ii+9,:)+3*sqrt(data.diag_P(9+ii,:)),'r:');
h(2)=plot(data.t, data.X(ii+9,:)-3*sqrt(data.diag_P(9+ii,:)),'r:');
legend(h,{'estimated','3 \sigma'})
ylabel('m/s^2')
xlabel('t [s]')
title('Acc bias')
grid minor;
end



figure(6)
clf
h=zeros(1,2);
for ii=1:3
subplot(3,1,ii)
hold on;
h(1)=plot(data.t, 180/pi*data.X(ii+12,:),'r');
h(2)=plot(data.t, 180/pi*(data.X(ii+12,:)+3*sqrt(data.diag_P(12+ii,:))),'r:');
h(2)=plot(data.t, 180/pi*(data.X(ii+12,:)-3*sqrt(data.diag_P(12+ii,:))),'r:');
legend(h,{'estimated','3 \sigma'})
ylabel('deg/s')
xlabel('t [s]')
title('Gyro bias')
grid minor;
end







end