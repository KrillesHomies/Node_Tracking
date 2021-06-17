function [T, zupt] = ZUPT(u,settings,detector_type)


% Allocate memmory
zupt=zeros(1,length(u));

% Run the desired detector type. Each detector return a vector with their 
% calculated test statistics T. 
switch detector_type
    
    case 'GLRT'
        T=GLRT(u,settings); 
    
    case 'MV'
        T=MV(u,settings);
        
    case 'MAG'
        T=MAG(u,settings);
        
    case 'ARE'
        T=ARE(u,settings);
        
    case 'PMV'
        T=PMV(u,settings);
        
    otherwise
        disp('The choosen detector type not recognized. The GLRT detector is used')
        T=GLRT(u,settings);
end

% Check if the test statistics T are below the detector threshold. If so, 
% chose the hypothesis that the system has zero velocity 
W=settings.window_size;
for k=1:length(T)
    if T(k)<settings.gamma
       zupt(k:k+W-1)=ones(1,W); 
    end    
end

% Fix the edges of the detector statistics
T=[max(T)*ones(1,floor(W/2)) T max(T)*ones(1,floor(W/2))];
end





%% SUBFUNCTIONS 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = GLRT(u)
%
%> @brief Function that runs the generalized likelihood test (SHOE detector). 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=GLRT(u,settings)


g=norm(settings.g);
sigma_a=settings.sigma_acc^2;
sigma_g=settings.sigma_gyro^2;
W=settings.window_size;


N=length(u);
T=zeros(1,N-W+1);

for k=1:N-W+1
   
    ya_m=mean(u(1:3,k:k+W-1),2);
    
    for l=k:k+W-1
        tmp=u(1:3,l)-g*ya_m/norm(ya_m);
        T(k)=T(k)+u(4:6,l)'*u(4:6,l)/sigma_g+tmp'*tmp/sigma_a;    
    end    
end

T=T./W;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = MV(u)
%
%> @brief Function that runs the acceleration moving variance detector. 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=MV(u,settings)



sigma_a=settings.sigma_acc^2;
W=settings.window_size;

N=length(u);
T=zeros(1,N-W+1);



for k=1:N-W+1
    
    ya_m=mean(u(1:3,k:k+W-1),2);
    
    for l=k:k+W-1
        tmp=u(1:3,l)-ya_m;
        T(k)=T(k)+tmp'*tmp;    
    end    
end

T=T./(sigma_a*W);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = MAG(u)
%
%> @brief Function that runs the acceleration magnitude detector. 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=MAG(u,settings)

g=norm(settings.g);
sigma_a=settings.sigma_acc^2;
W=settings.window_size;

N=length(u);
T=zeros(1,N-W+1);


for k=1:N-W+1
    for l=k:k+W-1
        T(k)=T(k)+(norm(u(1:3,l))-g)^2;    
    end    
end

T=T./(sigma_a*W);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = ARE(u)
%
%> @brief Function that runs the angular rate energy detector. 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=ARE(u,settings)



sigma_g=settings.sigma_gyro^2;
W=settings.window_size;

N=length(u);
T=zeros(1,N-W+1);


for k=1:N-W+1
    for l=k:k+W-1
        T(k)=T(k)+norm(u(4:6,l))^2;    
    end    
end

T=T./(sigma_g*W);

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = PMV(u)
%
%> @brief Function that cacluates the variance over the presure 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=PMV(u,settings)



sigma_p=settings.sigma_pressure^2;
W=settings.window_size;

N=length(u);
T=zeros(1,N-W+1);

for k=1:N-W+1
    
    ya_m=mean(u(1,k:k+W-1),2);
    
    for l=k:k+W-1
        tmp=u(1,l)-ya_m;
        T(k)=T(k)+tmp'*tmp;    
    end    
end

T=T./(sigma_p*W);

    
end
