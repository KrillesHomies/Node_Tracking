function [] = plot_trajectory(data, trajectory,NameGraph, figureNr,Color)

    figure(figureNr);
    sgtitle(NameGraph)
    % Plot
    subplot(3,1,1)
    plot(trajectory.GNSS.t(:),trajectory.GNSS.pos(1,:),'.k')
    title('North')
    hold on
    plot(data.t,data.X(1,:),'color',Color)
    ylabel('m')
    xlabel('t(s)')
    subplot(3,1,2)
    plot(trajectory.GNSS.t(:),trajectory.GNSS.pos(2,:),'.k')
    title('East')
    hold on
    plot(data.t,data.X(2,:),'color',Color)
    ylabel('m')
    xlabel('t(s)')
    subplot(3,1,3)
    plot(trajectory.GNSS.t(:),trajectory.GNSS.pos(3,:),'.k')
    title('Down')
    hold on
    plot(data.t,data.X(3,:),'color',Color)
    ylabel('m')
    xlabel('t(s)')
    
end