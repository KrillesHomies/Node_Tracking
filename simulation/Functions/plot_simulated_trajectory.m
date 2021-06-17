function [] = plot_simulated_trajectory(data)

    t = data.Data_N.IMU.t;
    figure(1);
    sgtitle('Position')
    % Plot
    subplot(3,1,1)
    plot(t,data.X_r(1,2:end),'.k')
    title('North')
    ylabel('m')
    xlabel('t(s)')
    subplot(3,1,2)
    plot(t,data.X_r(2,2:end),'.k')
    title('East')
    ylabel('m')
    xlabel('t(s)')
    subplot(3,1,3)
    plot(t,data.X_r(3,2:end),'.k')
    title('Down')
    ylabel('m')
    xlabel('t(s)')

    figure(2);
    sgtitle('Attitude')
    % Plot
    subplot(3,1,1)
    plot(t,data.X_r(7,2:end),'.k')
    title('Roll')
    ylabel('rad')
    xlabel('t(s)')
    subplot(3,1,2)
    plot(t,data.X_r(8,2:end),'.k')
    title('Pitch')
    ylabel('rad')
    xlabel('t(s)')
    subplot(3,1,3)
    plot(t,data.X_r(9,2:end),'.k')
    title('Yaw')
    ylabel('rad')
    xlabel('t(s)')
        
    
end