clc;
clear all;
close all;


n=5;
level=3;
display_inversion_circle = 1;
create_mpost = 0;
create_svg = 0;
filename_export = 'a1';

S = apollonian_2D(n,level,display_inversion_circle,...
                                create_mpost,create_svg,filename_export)