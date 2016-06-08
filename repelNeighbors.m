function [Dx, Dy] = repelNeighbors(idx, i, Cx, Cy)


tempX = Cx(idx) - Cx(i);
tempY = Cy(idx) - Cy(i);
dx = tempX./abs(tempX);
dy = tempY./abs(tempY);
% correct for case division by zero 
dx(isnan(dx)) = 0;
dy(isnan(dy)) = 0;

% % update contributions
Dx = -dx; 
Dy = -dy;


        
