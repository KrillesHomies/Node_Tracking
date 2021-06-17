%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [F G]=state_matrix(q,u)
%
%> @brief Function for calculating the state transition matrix F and 
%> the process noise gain matrix G. 
%>
%> @details Function for calculating the state transition matrix F and 
%> the process noise gain matrix G, given the current orientation of 
%> the platform and the specific force vector.  
%>
%> @param[out]   F     State transition matrix.
%> @param[out]   G     Process noise gain matrix.
%> @param[in]    u     IMU data [specific force, angular rates].
%> @param[in]    q     Old quaternions
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F G]=state_matrix(quat,u,dt)


R_B2N = q2dcm(quat);

St = skew(R_B2N*(u(1:3)'));

%[~,St] =Lineariztion_MAG_2(dcm2euler(q2dcm(quat)), u(1:3)');

O=zeros(3);

% Identity matrix
I=eye(3);

% Bias errors included
Fc=[O  I  O   O      O;
    O  O  St  R_B2N  O;
    O  O  O   O      -R_B2N;
    O  O  O   O       O;
    O  O  O   O       O];

% Noise gain matrix
Gc=[O O O O; R_B2N O O O; O -R_B2N O O; O O I O; O O O I];

% Approximation of the discret time state transition matrices
F=eye(15) + dt*Fc;
G= dt * Gc;
% G= dt * Gc;
end
