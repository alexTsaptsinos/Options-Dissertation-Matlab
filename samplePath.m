%% SAMPLE PATH GENERATOR
% This code generates a sample path for a 1D BM with the below parameters

steps = 40;
timestep = 1/steps;
W = 0;
points = [W];
mean = 0;
variance = 1;

for i = 1:steps;
    Z = randn;
    W = W + mean*timestep + sqrt(variance*timestep)*Z;
    points = [points W];
end

plot([0:steps],points);
  