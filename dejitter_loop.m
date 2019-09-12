clc;
clear all;
close all;
warning('off','all')
global kkk R2 tv_img_update
method = 3;     % 1 for fminunc, 2 for fminsearch, 3 for SA, 4 for ga integer
P = phantom('Modified Shepp-Logan',200);
figure(1)
imshow(P)
kk=1;
ths_end1 = 90;
ths_cut = 4;

theta = 0:kk:ths_end1;
n_s = ths_cut/kk;
[R0,xp] = radon(P,theta);
k_loop = 1;
for i=1:ths_cut:ths_end1-ths_cut
r= randi([-20 20],1,n_s+1);
k_inner = 1;
for ii = 1:n_s+1
    R1(:,ii) = imtranslate(R0(:,i+k_inner),[0,r(ii)]);
    k_inner = k_inner + 1;
end
figure(3)
imagesc(R1)
kkk=1;
x0 = zeros(1,n_s+1);
if (method ==1)
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',30000,'MaxFunctionEvaluations',90000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
problem = createOptimProblem('fmincon','x0',x0, ...
    'objective',@(x)shiftingpixel(x,R1,kkk),'options',opts);
[x_soln] = fmincon(problem);
elseif (method ==2)   
opts = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',10^(-10),'TolX',10^(-10));
[x_soln] = fminsearch(@(x)shiftingpixel(x,R1,kkk),x0,opts);
elseif (method ==3)
    opts = optimoptions('simulannealbnd','MaxIter',3000,'PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplottemperature});
    lb = -50*ones(1,n_s+1);
ub = -lb;
    x_soln = simulannealbnd(@(x)shiftingpixel(x,R1,kkk),x0,lb,ub,opts);
elseif (method ==4)
opts = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
x_soln = ga(@(x)shiftingpixel(x,R1,kkk),n_s+1,[],[],[],[],[],[],[],opts);
end
R_assemble(:,:,k_loop) = R2;
k_loop = k_loop + 1;
end

figure(10)
aa=size(R_assemble);
R_final = reshape(R_assemble,aa(1),aa(2)*aa(3));
imagesc(R_final)

