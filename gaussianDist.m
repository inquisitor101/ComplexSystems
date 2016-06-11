function out = gaussianDist(SD, N)

% SD: standard deviation
% mu: mean/average

mu = 0;
out = normrnd(mu, SD, [N, 1]);

