clc;
clear all;
choice = 'dist';    % man for manuel, gmsh for gmsh
filename1 = 'gmsh_in1.txt';
filename2 = 'gmsh_in2.txt';

if (strcmp(choice,'man'))
n2 = 10;             % Internel mesh number angle wise
n3 = 10;            % Internel mesh number radius wise

hh = r/n3;
h2 = 2*pi/n2;
rho = 0:hh:r;
theta1 = -pi : h2 : pi;
[th_msh, rr_msh] = meshgrid(theta1, rho(1:end-1));
x_fxloc = rr_msh.*cos(th_msh);
y_fxloc = rr_msh.*sin(th_msh);

figure(3)
s1 = size(x_fxloc,1);
s2 = size(x_fxloc,2);
sll = s1 +ss;
x_plot = reshape(x_fxloc,[],1);
y_plot = reshape(y_fxloc,[],1);
scatter(x_plot,y_plot);
save('geo2.mat')
elseif (strcmp(choice,'gmsh'))
    mesh = load(filename1);

    x_plot = mesh(:,2);
    y_plot = mesh(:,3);
%     [X_img,Y_img] = meshgrid(x_plot,y_plot);
    tri = delaunay(x_plot,y_plot);

%     scatter(mesh(:,1),mesh(:,2),20,rand(66,1),'filled');
 save('geo2.mat')
 
elseif  (strcmp(choice,'dist'))
    fd=@(p) sqrt(sum(p.^2,2))-1;
    [p_plot,tri]=distmesh2d(fd,@huniform,0.1,[-1,-1;1,1],[]);
    x_plot = p_plot(:,1);
    y_plot = p_plot(:,2);
    scatter(x_plot,y_plot);
     save('geo2.mat')
end


%%%% Set data measurement location:
