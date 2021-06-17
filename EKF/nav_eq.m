%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navigation equations using the 
%
% Input -   X - State (Position, Velocity, Euler Angles)[1x9]
%           u_h - Accelerometer and Gyroscope Data [1x6]
%           quat - Orientation expressed in quaternions
%           g - gravity vector
%           dt - timestep
%
% Ouput - X - Updated State 
%         quat - Updated orientation in quaternoins
%
% Edited by Kristoffer Lindve 2021-06-08
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,quat] = navigation_eq(X,u_h,quat,g,dt)

    % Motion model
    a_meas = u_h(1:3); % measure Acceleration
    ang_measured = u_h(4:6); % Angualr Velocity
    
    X(1:6)=update_linear_movment(X(1:6),a_meas,quat,g,dt); % Apply linear movment
    quat=update_orientation(ang_measured,quat,dt); %Get new Rotation
    
    %Update state
%     R_N2B = R_B2N';
    X(7:9) = dcm2euler(q2dcm(quat)); %X(7:9) = rotm2eul(R_N2B,'XYZ');
    


function X=update_linear_movment(X,a_meas,quat,g,dt)

A = [1 0 0 dt 0 0; ...
     0 1 0 0 dt 0; ...
     0 0 1 0 0 dt; ...
     0 0 0 1 0 0; ...
     0 0 0 0 1 0; ...
     0 0 0 0 0 1];
B = [0.5*dt^2 0        0; ...
     0        0.5*dt^2 0; ...
     0        0        0.5*dt^2; ...
     dt       0        0; ...
     0        dt       0; ...
     0        0        dt];

a_N = q2dcm(quat)*(a_meas) - g'; % Calucalte Acceleration

X(1:6) = A*X(1:6) + B*a_N; % Update Position anc Velocity

return


function quat=update_orientation(ang_vel_meas,quat,dt)

    ang_vel = ang_vel_meas;
    P=ang_vel(1)*dt;
    Q=ang_vel(2)*dt;
    R=ang_vel(3)*dt;

    OMEGA=zeros(4);
    OMEGA(1,1:4)=0.5*[0 R -Q P];
    OMEGA(2,1:4)=0.5*[-R 0 P Q];
    OMEGA(3,1:4)=0.5*[Q -P 0 R];
    OMEGA(4,1:4)=0.5*[-P -Q -R 0];

    v=norm(ang_vel)*dt;
    if v~=0
        quat=(cos(v/2)*eye(4)+2/v*sin(v/2)*OMEGA )*quat;
        quat=quat./norm(quat);
    end
    
return