clc;
clear all;

data = hdf5read('i23_12641.h5','/data_raw');
figure(1)
A= squeeze(data(:,400,:));   
imagesc(A);

figure(2)
B= squeeze(data(100,:,:));
imagesc(B);
% % T = randn()