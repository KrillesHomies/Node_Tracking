%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman Smoother smooths the data and outputs the smoothed state. 
%
% Input -   indata - Output from the kalman filter
%
% Ouput - outdata - Updated path with updated process covariance 
%
% Edited by Kristoffer Lindve 2021-06-08
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [outdata] = RTS(indata,settings)

itterations = size(indata.X,2);

X = indata.X;
P = indata.P;
x = indata.X_predicted;
p = indata.P_predicted;
F = indata.F;
%H_mag = [zeros(3,6) skew(q2dcm(x_h(7:10))*settings.MagneticField') zeros(3,6)];
for i= itterations-1:-1:1
        
    % Get F and C matrix 
    quat = dcm2q(Rt2b(X(7:9,i))');
    %[F G]=state_matrix(quat,indata.u_h(:,i)',settings.dt);

    C = P(:,:,i)*F(:,:,i+1)'*pinv(p(:,:,i+1));
    %calcualte innovation
    y = [[X(:,i+1)] - [x(:,i+1)]];
    z = C*y;
    X(:,i) = X(:,i) + z;
    % Update Angles
%     quat=Gamma(quat,z(7:9)');
%     X(7:9,i) = dcm2euler(q2dcm(quat/norm(quat)));
   
    % Update Covariance
    P(:,:,i) =  P(:,:,i) + C*(P(:,:,i+1) - p(:,:,i+1))*C';
    P(:,:,i) = (P(:,:,i) + P(:,:,i)')/2; 

end

outdata = indata;
outdata.X = X;
outdata.P = P;

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