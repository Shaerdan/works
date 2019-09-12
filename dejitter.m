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
ths_end1 = 180;
ths_cut = 10;
theta = 0:kk:ths_end1;
n_s = ths_cut/kk;
[R0,xp] = radon(P,theta);

% Igrad0a = imgradient(R0);
[gx0,gy0] = imgradientxy(R0);
tv_img0a = sum(sum((sqrt(gx0.^2+gy0.^2))));
% tv_img0a = norm(Igrad0a,2);

R = R0(:,80:ths_cut+80);
Igrad0b = imgradient(R);
[gx1,gy1] = imgradientxy(R);
tv_img0b = sum(sum((sqrt(gx1.^2+gy1.^2))));

figure(2)
imagesc(R)
% r = 5*randn(1,n_s+1);
r= randi([-20 20],1,n_s+1);

for i = 1:n_s
    R1(:,i) = imtranslate(R(:,i),[0,r(i)]);
end

R0 = R1;

figure(3)
imagesc(R1)

kkk=1;
x0 = zeros(1,n_s+1);
% lb = -ones(1,n_s+1);
% ub = -lb;
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
    opts = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplottemperature});
    
%     opts = optimset('Display','iter','PlotFcns',@optimplotfval);

lb = -50*ones(1,n_s+1);
ub = -lb;
    x_soln = simulannealbnd(@(x)shiftingpixel(x,R1,kkk),x0,lb,ub,opts);
elseif (method ==4)
opts = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
x_soln = ga(@(x)shiftingpixel(x,R1,kkk),n_s+1,[],[],[],[],[],[],[],opts);

end

figure(4)
imagesc(R1)

figure(5)
imagesc(R1)

figure(6)
imagesc(R2)

