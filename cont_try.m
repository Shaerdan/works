close all;
[xq,yq] = meshgrid(-2:.005:2, -2:.005:2);
vq = griddata(mesh(:,1),mesh(:,2),R',xq,yq);

contour(xq,yq,vq,30)
hold on
% plot3(mesh(:,1),mesh(:,2),R','*')
% xlim([-2.7 2.7])
% ylim([-2.7 2.7])