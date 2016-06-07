function dist = getDistance(i, L, Cx, Cy, N)


A  = zeros(N, 9);   % 9 different position scenarios
% difference between (N-1) and current position
xR = Cx(i) - Cx(:);
yR = Cy(i) - Cy(:);
% position scenarios
A(:,1) = sqrt((xR)  .^2 + (yR)  .^2);
A(:,2) = sqrt((xR-L).^2 + (yR)  .^2);
A(:,3) = sqrt((xR)  .^2 + (yR-L).^2);
A(:,4) = sqrt((xR+L).^2 + (yR)  .^2);
A(:,5) = sqrt((xR)  .^2 + (yR+L).^2);
A(:,6) = sqrt((xR+L).^2 + (yR+L).^2);
A(:,7) = sqrt((xR+L).^2 + (yR-L).^2);
A(:,8) = sqrt((xR-L).^2 + (yR+L).^2);
A(:,9) = sqrt((xR-L).^2 + (yR-L).^2);

A(i, :) = L; % L is max, so it will never affect solution

dist = A';