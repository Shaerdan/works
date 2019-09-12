clc;
clear all;
close all;
%%%%% ERROR CONTOUR R = l2norm(u_bc(xs,ys) - u_bc(xi,yi)): 
global flag1 flag2

load('geo.mat')
load('data.mat');
sanity = 0;
%  figure(6)    
% Trigular mesh: 
% load('contmesh.mat');
% x_source0(:,1) = p(:,1);
% x_source0(:,2) = p(:,2);

% Structured mesh:
load('mesh.mat')
x_source0(:,1) = mesh(:,1);
x_source0(:,2) = mesh(:,2);


phim_noisy = phim_measure_noisy;

m = length(x_source0);
k = 1;
for i = 1:dim:(m*dim)
x_loc(i) = x_source0(k,1);
x_loc(i+1) = x_source0(k,2);
k = k + 1;
end

    k = 1; 
    s = 1;
for ii = 1:dim:(m*dim) 
     [ phim,pkm,F1 ] = Sol_Poisson_BEM( x_loc(ii:ii+1),S,x,y,xx,yy,nx,ny,xm,ym,lx,...
                                   ly,l,r,theta,sanity,dim);
                        
    for m=1:length(xd)
          F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
          phim_m(m) =  F(2,:)*phim' - F(1,:)*pkm';
    end
    
%     norm(phim_m,2);
    R(k) = norm((phim_m-phim_noisy),2);
        if ( k > 1 ) 
            dR = R(k)/R(k-1);
            if (dR > 100)
            mark_x(s) = x_loc(ii);
            mark_y(s) = x_loc(ii+1);
            val(s) = dR;
            s = s + 1;
            end
        end

    k = k + 1 ;
end
%             save('debug.mat','mark_x','mark_y','val')

% z = linsolve(A,b);
% Other linear solvers:

% R(k:m) = R(k);
px = 0.001;
[xq,yq] = meshgrid(-1:px:1, -1:px:1);
vq = griddata(mesh(:,1),mesh(:,2),R',xq,yq);
figure(7)
contour(xq,yq,vq,50); hold on;

load('inv_soln.mat')
plot(x,y); hold on;
% contourf(mesh(:,1),mesh(:,2),R)
ll = length(x_source);
for i = 1:dim:ll
scatter(x_source(i),x_source(i+1),'b','*'); hold on; 
scatter(loc_inv(i),loc_inv(i+1),'b','+'); hold on;
end


% contourf(t,p(:,1),p(:,2),R','FaceColor','interp','EdgeColor','interp')
% [cout,hout] = tricontour(p,t,R',1000);
% hold on; 
% scatter(xcc,ycc)
