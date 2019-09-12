clc;
clear all;
fd=@(p) sqrt(sum(p.^2,2))-0.998;
[p,t]=distmesh2d(fd,@huniform,0.025,1.1*[-0.6,-0.6;0.6,0.6],[]);
save('contmesh.mat')
1331