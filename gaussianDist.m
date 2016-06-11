function out = gaussianDist(SD, N)

% SD: standard deviation
% mu: mean/average

mu = 0;
out = SD.*randn(N,1) + mu;

