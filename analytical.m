clc;
clear all;
close all;

swtch = 1;
if (swtch == 1)
load('data.mat');
end
if (swtch == 2)
load('datafixloc.mat')
end
if (swtch == 3)
load('datafixloc.mat')    
    load('gain.mat')
end

n_th = 2000;
ep = 0.0001;
h = (2*pi)/(n_th-1);
theta = 0:h:2*pi;
sanity = 0;
n_n = 1000;
global alph;
alph = 2;
%%%% commpute vector g:
    for i=1:n_source
    temp1 = (cos(theta)-x_plot(i)).^2 + (sin(theta)-y_plot(i)).^2;
%     pk(k,:) = (1/(2*pi))-(S(k)/(2*pi))*gx;
    g_th(i,:) = - (S(i)/(2*pi))*(1./(sqrt(temp1)));
    end
    g_th =  sum(g_th,1);
%     figure(1)
%     plot(g);
%     figure(2)
%     plot(f_data);

%%%%%% SANITY TEST %%%%%%
if (sanity == 1)
g_th = 2*((cos(theta)).^2 - (sin(theta)).^2);
end
if (sanity == 2)
   g_th= (alph*1^(alph-1))*cos((alph+1)*theta);
end
 
    
for i = 1:n_th
    for j = 1:n_n
        CS(i,j) = j*cos((j)*theta(i));
        SI(i,j) = j*sin((j)*theta(i));
    end
end
  
K = [CS SI];

s = svd(K);
figure(1)
plot(log(s));
lam = 0.0001;
coef0 = pinv(K)*g_th';
coef1 = (K'*K+lam*eye(2*n_n,2*n_n))\(K'*g_th');

A1 = coef1(1:floor(end/2));
B1 = coef1(floor(end/2)+1:end);

A0 = coef0(1:floor(end/2));
B0 = coef0(floor(end/2)+1:end);

err = (coef0 - coef1);

figure(2)
plot(A1); hold on; plot(A0);

figure(3)
plot(B1); hold on; plot(B0);
save('analytical.mat','A1','B1','A0','B0');

r = 1; 
n_mesh = 100;
h_mesh = 2*pi/n_mesh;
theta_mesh = (0:h_mesh:2*pi);
hh = r/10;
rho = 0:hh:r;
[th_mesh, rr_mesh] = meshgrid(theta_mesh, rho);
xc = rr_mesh.*cos(th_mesh);
yc = rr_mesh.*sin(th_mesh);

%  for i=1:size(xc,2)-1
%      xcm(:,i) = 0.5*(xc(:,i) + xc(:,i+1));
%      ycm(:,i) = 0.5*(yc(:,i) + yc(:,i+1));
%  end


for n = 1:length(A1)
   ser1 = A1(n)*cos(n*th_mesh);
   ser2 = B1(n)*sin(n*th_mesh);
   ser = ser1 + ser2;
   u_anl(:,:,n) = (rr_mesh.^n).*ser;
end
    
u_anl = sum(u_anl,3);
save('analytical.mat','A1','B1','A0','B0','u_anl');
figure(4)
surf(rr_mesh,th_mesh,u_anl);
%        zlim([-1 1])


[ u_bem ] = solve_laplace_neumann_F( th_mesh, rr_mesh,xc,yc,sanity,swtch );
if (sanity == 2)
    u_bem = (rr_mesh.*cos(th_mesh)).^2 -(rr_mesh.*sin(th_mesh)).^2 ;
end
figure(5)
surf(rr_mesh,th_mesh,u_bem);
%        zlim([-1 1])
       
       figure(6)
       subplot(2,2,1)
       surf(xc,yc,u_anl)
       subplot(2,2,2)
       surf(xc,yc,u_bem)


       
