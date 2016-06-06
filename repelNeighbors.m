function [Dx, Dy] = repelNeighbors(N, i, Cx, Cy)

dx = 0; dy = 0;  % reset
for j=1:i-1  % unroll loop 
    tempX = Cx(j) - Cx(i); % difference X 
    tempY = Cy(j) - Cy(j); % difference Y
    dx = dx + tempX/abs(tempX); % normalized X
    dy = dy + tempY/abs(tempY); % normalized Y
end
for j=i+1:N  % unroll loop 
    tempX = Cx(j) - Cx(i); % difference X
    tempY = Cy(j) - Cy(j); % difference Y 
    dx = dx + tempX/abs(tempX); % normalized X
    dy = dy + tempY/abs(tempY); % normalized Y
end
% update contributions
Dx = -dx; 
Dy = -dy;
        
