clc;
clear all;
close all;

global k_star ak bk n_source lambda S_real ths_real rs_real
lambda = 10^(-6);

r=0.6;
% ths = 0.01;
I = 1;
k_star = 1000;

samples = k_star;
n_spec = 30;
level = (10^-2);
n_fft = samples;
n_sample = samples;
h_integral = 2*pi/(n_sample-1);
kk = 1;

    theta = -pi:h_integral:pi;

    S_real = [1 -1 ];
    n_source = length(S_real);
    ang1 = 0.846; small_ang = 0.02;
    ang2 = ang1 + small_ang;
    fan_range = ang2 - ang1;
    ang_stp = (fan_range)/(n_source-1);
    ths_real = ang1:ang_stp:ang2;
    rs_real = r;

rs_real = rs_real*ones(1,n_source);

for i = 1:n_source 
phi_real(:,i) = (S_real(i)/(2*pi))*log(sqrt(1+rs_real(i)^2 ...
    -2*rs_real(i)*cos(theta-ths_real(i))));
end

phi_real = 2*sum(phi_real,2);


for i = 1:k_star
    func1 = @(xx) cos(i*xx).*phi_func(xx);
    func2 = @(xx) sin(i*xx).*phi_func(xx);
    an2(i) = (1/pi)*integral(func1,-pi,pi);
    bn2(i) = (1/pi)*integral(func2,-pi,pi);
end

    for i = 1:k_star
        ak1(i) = -(S_real(1)/pi)*((rs_real(1)^i)/i)*cos(i*(ths_real(1)));
        ak2(i) = -(S_real(2)/pi)*((rs_real(2)^i)/i)*cos(i*(ths_real(2)));
        bk1(i) = -(S_real(1)/pi)*((rs_real(1)^i)/i)*sin(i*(ths_real(1)));
        bk2(i) = -(S_real(2)/pi)*((rs_real(2)^i)/i)*sin(i*(ths_real(2)));
    end
    
    ak = ak1+ak2;
    bk = bk1+bk2;


% ak = an2(2:k_star); bk = bn2(2:k_star);

%%% 1. nutral source strength constraint 
temp_A = [1 0 0];
Aa = [repmat(temp_A,1,n_source)];    
ba = [0];

%%% 2. upper and lower bound for source strength constraint
K=2*max(S_real); temp_A1 = [1 1/K pi/K]; temp_A2 = [1 1/K pi/K];

lb = [-K*repmat(temp_A1,1,n_source)]'; ub =-lb;
[K*repmat(temp_A2,1,n_source)]';


x0=-4:(8/(3*n_source-1)):4;

opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',30000,'MaxFunctionEvaluations',90000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16);
problem = createOptimProblem('fmincon','x0',x0, ...
    'objective',@fourier_obj_fun,'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,...
    'options',opts);



x_soln = fmincon(problem);

k=1;
for i=1:n_source
S_soln(i) = x_soln(k);
rs_soln(i) = x_soln(k+1);
ths_soln(i) = x_soln(k+2);
k=k+3;
end
 xs_soln = rs_soln.*cos(ths_soln); ys_soln = rs_soln.*sin(ths_soln);

figure(11)
 xs = rs_real.*cos(ths_real); ys = rs_real.*sin(ths_real);
 x = cos(theta); y1 = sin(theta);
 plot(x,y1); hold on;
    k = 1;
    for i = 1:n_source
        scatter(xs(i),ys(i));
        text1 = S_real(k);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        k = k + 1;
        text_fin =text(xs(i)+0.05,ys(i)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 12;
        hold on;
    end

figure(12)
 plot(x,y1); hold on;
    k = 1;
    for i = 1:n_source
        scatter(xs(i),ys(i),'*');
        scatter(xs_soln(i),ys_soln(i));
        text1 = S_soln(k);
        text1 = num2str(text1);
        text1 = cellstr(text1);
        k = k + 1;
        text_fin =text(xs_soln(i)+0.05,ys_soln(i)+0.05,text1);
        fontsize = text_fin.FontSize;
        text_fin.FontSize = 10;
        hold on;
    end


    
%     for i = 1:k_star-1
%         cosk1(i,:) = cos(i*theta);
%         sink1(i,:) = sin(i*theta);
%     end
%     ak1 = ak1';
%     ak2 = ak2';
% 
%     bk1 = bk1';
%     bk2 = bk2';
%     cosk1 = cosk1';    
%     sink1 = sink1';
%     
%     for i = 1:length(theta)
%     phi_fourier(i) = sum(ak1.*cosk1(i,:)') + sum(ak2.*cosk1(i,:)')  +...
%          sum(bk1.*sink1(i,:)') + sum(bk2.*sink1(i,:)');
%     end
%     
%     for i = 1:length(an2)
%         coskd1(i,:) = cos(i*theta);
%         sinkd1(i,:) = sin(i*theta);
%     end
%     
%     coskd1 = coskd1';
%     sinkd1 = sinkd1';
% %     
%     for i = 1:length(theta)
%         phi_fourierd(i) = sum(an2.*coskd1(i,:)) + sum(bn2.*sinkd1(i,:));
%     end
% 
%     
%     
%     figure(5)
%     plot(theta,phi_fourier,'*'); hold on;  plot(theta,phi_real);
%     
%     figure(6)
%     plot(theta,phi_fourier,'*'); 
%     
%     figure(7)
%     plot(theta,phi_fourierd,'*'); 