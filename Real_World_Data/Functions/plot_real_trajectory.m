%% Plots the trajectory with the 


function [] = plot_real_trajectory(data,time_start, gps_data,drop, pick,NameGraph, figureNr)

     ellipse = @(point, n_unc, e_unc) [n_unc*cos(-0.1:0.1:(2*pi)); e_unc*sin(-0.1:0.1:(2*pi))]'+point;

     figure(figureNr);
     title(NameGraph)
     hold on
     %% GPS
     %scatter3(gps_data.GNSS.pos(1,:),gps_data.GNSS.pos(2,:),-gps_data.GNSS.pos(3,:),'.k')
     colors = [1 0 0; 0 1 0; 0 0 1];
     %% Drop point
     h(1) = scatter3(drop(1),drop(2),-drop(3),'MarkerFaceColor',colors(1,:));
     e = ellipse(drop(1:2),sqrt(drop(4)),sqrt(drop(5))); 
     plot3((e(:,1)),(e(:,2)),-drop(3)*ones(size(e(:,1))),'-.','color',colors(1,:));
     %% Pick up 
     h(2) = scatter3(pick(1),pick(2),-pick(3),'MarkerFaceColor',colors(2,:));
     e = ellipse(pick(1:2),sqrt(pick(4)),sqrt(pick(5))); 
     plot3((e(:,1)),(e(:,2)),-pick(3)*ones(size(e(:,1))),'-.','color',colors(2,:));
     %% Predicted Position 
     h(3) = scatter3(data.X(1,end),data.X(2,end),-data.X(3,end),'MarkerFaceColor',colors(3,:));
     e = ellipse(data.X(1:2,end)',sqrt(data.P(1,1,end)),sqrt(data.P(2,2,end))); 
     plot3((e(:,1)),(e(:,2)),-data.X(3,end)*ones(size(e(:,1))),'-.','color',colors(3,:));
     %% Smoothed Predicted trajectory
     [~,n] = min(abs(data.t - time_start));
     
     plot3(data.X(1,n:end),data.X(2,n:end),-data.X(3,n:end))
     
     %% Plot orientation
     L = size(data.X(1,n:end),2);
     nr_coord_systems = 100;
     coords = round(linspace(L,nr_coord_systems)+n-1,0);
     
     for i=1:nr_coord_systems
         
         RN2B = Rt2b(data.X(7:9,coords(i))');
         color = 'rgb';
         for j=1:3
             quiver3([data.X(1,coords(i))]',[data.X(2,coords(i))]',-[data.X(3,coords(i))]',RN2B(j,1),RN2B(j,2),-RN2B(j,3),color(j),'LineWidth',1);
         end
        
     end
     
     %% Labels
     xlabel('North (m)')
     ylabel('East (m)')
     zlabel('Up (m)')
     
     legend(h, {'Drop-off','Pickup','Predicted'})
end