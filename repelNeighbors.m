function [Dx, Dy] = repelNeighbors(idx, i, Cx, Cy)


tempX = Cx(idx) - Cx(i);
tempY = Cy(idx) - Cy(i);
dx = tempX./sqrt(tempX.^2 + tempY.^2);
dy = tempY./sqrt(tempX.^2 + tempY.^2);
% correct for case division by zero 
dx(isnan(dx)) = 0;
dy(isnan(dy)) = 0;

% % update contributions
Dx = -sum(dx); 
Dy = -sum(dy);


        
