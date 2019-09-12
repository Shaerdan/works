clc;
clear all;
close all;


r=1;               %radius of the disk
k = 1;
% alph = 1;          %dimention parameter for periodic bc in sanity test
% levels = 40;       %error contour levels
% % Soure term parameters, take S = 0 for no source term case:
% % Switches:
% type = 'N';     % 1.Boundary type switch,  D--> Dirichlet, N --> Neumann.  
% sanity = 0;     % 2.Sanity Test Switch, if sanity == 1, perform u = u(r) = r^(alpha)*cos(alpha*theta)
% % test, Dirichlet boundary condition is u = u(r_bc), in Neumann bc is
% % du(r)/dr at r_bc;
% %%control pamameters for gmres solver:
% gmres_restart = 25;
% gmres_tol = 10^(-10); 
% gmres_maxit = 1000;

global alph snr flag1 
snr = 'no';          %signal to noise ratio for artificial noisy boundary%
alph = 2;          %dimention parameter for periodic bc in sanity test
% alph_tik = 0.00001;   % Tichnov Regularisation for NtD map
alph_deflate = 10^8; % deflation of NtD map

lambda_inv = 0.001; % Tichnov Regularisation for the inverse problem
data_pts = 100;
sanity = 0;
data_perccet=1/data_pts;      
full_partial = 0;   % 1 for full data, 0 for partial;
manuel = 1;  

       n_sp = 1;    % positive charges
       n_sm = 1; 
%        S1 = eye(n_sp,1);
%        S2 = -fliplr(S1);
        rr = 0.7;
        S = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];

%%%%%%%%%%% random angular location %%%%%%%%   
     fan_range = 2*pi;  % scattered range angle 
    n_source = length(S);
%     h = fan_range/n_source;
%        ang_stp = (2*pi)/n_source;
%         th = 0:ang_stp:2*pi;
     th = fan_range*rand(n_source,1);        

        for i = 1:n_source
         x_plot(i) = rr*cos(th(i));
         y_plot(i) = rr*sin(th(i));
        end 


msh = 10:100:1500;
  
for i = msh
r=1;               %radius of the disk
alph = 1;          %dimention parameter for periodic bc in sanity test
levels = 40;       %error contour levels
% Soure term parameters, take S = 0 for no source term case:
% Switches:
type = 'N';     % 1.Boundary type switch,  D--> Dirichlet, N --> Neumann.  
sanity = 0;     % 2.Sanity Test Switch, if sanity == 1, perform u = u(r) = r^(alpha)*cos(alpha*theta)
% test, Dirichlet boundary condition is u = u(r_bc), in Neumann bc is
% du(r)/dr at r_bc;
%%control pamameters for gmres solver:
gmres_restart = 500;
gmres_tol = 10^(-14); 
gmres_maxit = 10000;

global alph snr flag1 
snr = 'no';          %signal to noise ratio for artificial noisy boundary%
alph = 2;          %dimention parameter for periodic bc in sanity test
alph_tik = 0.000001;   % Tichnov Regularisation for NtD map
alph_deflate = 10^8; % deflation of NtD map

lambda_inv = 0.000001; % Tichnov Regularisation for the inverse problem
data_pts = 100;
sanity = 0;
data_perccet=1/data_pts;      
full_partial = 0;   % 1 for full data, 0 for partial;
manuel = 1;  

global alph 
alph = 1;          %dimention parameter for periodic bc in sanity test
levels = 40;       %error contour levels
% Soure term parameters, take S = 0 for no source term case:
% Switches:
sanity = 0;  
    
    
n1=i;             %boundary mesh number
h1=2*pi/n1;
theta=-pi:h1:pi;
[ x,y,xx,yy ] = mesh_circle( r,theta );
%check the plot
[ nx,ny,xm,ym,lx,ly,l ] = geom( x,y );



%%%% data generate %%%%%%



%         th = [0,ang_stp,2*ang_stp,3*ang_stp,3.1*ang_stp,4*ang_stp,2*pi];
       


    [ g, gk ] = boundary_pk( r,theta,n_source,nx,ny,x,y,S,x_plot,y_plot,sanity );

    [ ff,F1,A,P ] = NtD( g,xx,yy,...
                              nx,ny,xm,ym,lx,...
                              ly,l,alph_tik);

%     stp = data_perccet*n1;
%     xd = xm(1:stp:end);
%     yd = ym(1:stp:end);
%     
%     [ ve ] = get_ve( xd,yd,S, x_plot, y_plot );
% %     for m=1:length(xm)                               
%     for m=1:length(xd)
%           F = integrals(l,nx,ny,xx,yy,xd(m),yd(m),lx,ly);
%           f_data(m) =  2*(F(2,:)*ff' - F(1,:)*g');
%           M(m,:) = F(2,:);
%           N(m,:) = F(1,:);
%     end  
% 
%     f_data = f_data + ve;
    
  %%%%%%%%% on boundary %%%%%%%%%%  
    for m=1:length(xm)
          F2 = integrals(l,nx,ny,xx,yy,xm(m),ym(m),lx,ly);
          f_data_reconst(m) =  2*(F2(2,:)*ff' - F2(1,:)*g');
          M(m,:) = F1(2,:);
          N(m,:) = F1(1,:);
    end  
    %%%%%% compute analytical solution %%%%%%
    for i = 1:n_source
    phi(:,i) =(S(i)/(2*pi))*log(sqrt((xm-x_plot(i)).^2 +...
        (ym-y_plot(i)).^2));
    end
    phi = sum(phi,2);
    error_cvg = (f_data_reconst - phi');
    record(k) = norm(error_cvg,2);
    record2(k) = norm(error_cvg,2)/(sqrt(length(N)));
    k = k + 1 ;
    
%     
%        [ ve ] = get_ve_forward( xc,yc,S,x_plot,y_plot );
% 
%        soln_poisson = u + ve;   
    clearvars -except record record2 r k n_source th x_plot y_plot S msh
   
end



n1=100;             %boundary mesh number
h1=2*pi/n1;
theta=-pi:h1:pi;
[ x,y,xx,yy ] = mesh_circle( r,theta );
%check the plot
[ nx,ny,xm,ym,lx,ly,l ] = geom( x,y );
            figure(1)
        plot(x,y); hold on;
    k = 1;
    for i = 1:n_source
        scatter(x_plot(i),y_plot(i));
        text1 = S(k);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        k = k + 1;
        text_fin =text(x_plot(i)+0.05,y_plot(i)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 4;
        hold on;
    end

     figure(2)
       loglog(msh,record2);
       xlabel('number of mesh points');
       ylabel('RMS error');