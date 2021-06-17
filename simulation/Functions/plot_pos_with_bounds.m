function [] = plot_pos_with_bounds(data, true_path,point,NameGraph, figureNr,Color)

ellipse = @(point, n_unc, e_unc) [n_unc*cos(-0.1:0.1:(2*pi)); e_unc*sin(-0.1:0.1:(2*pi))]'+point;

figure(figureNr)
scatter(true_path.X_r(1,point),true_path.X_r(2,point),500,'*k','LineWidth',3);
hold on 
scatter(data.X(1,point),data.X(2,point),'MarkerFaceColor',Color);
e = ellipse(data.X(1:2,point)',sqrt(data.P(1,1,point)),sqrt(data.P(2,2,point))); 
plot((e(:,1)),(e(:,2)),'-.','color',Color);
xlabel('North')
ylabel('East')
title(NameGraph)
    
    
end