function [  ] = plotsoln( x_invsoln,n_source,dim,x,y,x_source )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
S_inv = x_invsoln(1:n_source);
loc_inv = x_invsoln(1+n_source:end);
plot(x,y);
hold on;
for i = 1:dim:length(loc_inv)
scatter(x_source(i),x_source(i+1),'r','o'); 
hold on; 
end

for i = 1:dim:length(loc_inv)
scatter(loc_inv(i),loc_inv(i+1),'k','*')
hold on;
end

end

