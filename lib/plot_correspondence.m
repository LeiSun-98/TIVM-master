function plot_correspondence(pts_3d,pts_3d_,outlier_ratio)


figure;

n_ele=size(pts_3d,1);

pc1=pointCloud(pts_3d(1:n_ele,:));
pc2=pointCloud(pts_3d_(1:n_ele,:));

pcshow([pc1.Location(:,1),pc1.Location(:,2),pc1.Location(:,3)],[0 0 1],'MarkerSize',90);

hold on;

pcshow([pc2.Location(:,1),pc2.Location(:,2),pc2.Location(:,3)],[1 0 1],'MarkerSize',200);

for i=round(n_ele*outlier_ratio)+1:n_ele
    
    plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'g','LineWidth',1);
    
end

for i=1:round(n_ele*outlier_ratio)
    
    pp=plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'r','LineWidth',1);

    pp.Color(4) = 0.6;
    
end

grid off;
axis off;

end

