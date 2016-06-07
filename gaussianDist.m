function out = gaussianDist(range, stdDev, N)

out = zeros(N, 1); % make sure its an N-by-1
% step 1: generate numbers
a = range(1); 
b = range(2);
x = a + (b-a)*rand(N, 1);
mu = (a+b)/2; 


% step 2: gaussian function
coef1 = 1.0/( stdDev*sqrt(2*pi) );
coef2 = ( (x - mu)/stdDev );
coef2 = coef2.*coef2;

out = coef1.*exp(-coef2/2);
