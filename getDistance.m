function dist = getDistance(i, L, Cx, Cy, N)

% dist = zeros(N-1, 1); % initialize distance
% 
% for j=1:i-1
%     x = Cx(i) - Cx(j);
%     y = Cy(i) - Cy(j);
%     dist(j) = sqrt( x*x + y*y );
% end
% 
% for j=i+1:N
%     x = Cx(i) - Cx(j);
%     y = Cy(i) - Cy(j);
%     dist(j-1) = sqrt( x*x + y*y );
% end


% % % % % % % % % % % % % % % % % % % % % % % % %

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

A(i, :) = [];

dist = A';