function [eul] = dcm2euler(R)

% Get the corrected roll, pitch and heading from the corrected rotation
% matrix
eul(1)=atan2(R(3,2),R(3,3));
eul(2)=-atan(R(3,1)/sqrt(1-R(3,1)^2));
eul(3)=atan2(R(2,1),R(1,1));
