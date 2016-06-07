function [Dx, Dy] = repelNeighbors(idx, i, Cx, Cy)


tempX = Cx(idx) - Cx(i);
tempY = Cy(idx) - Cy(i);
dx = tempX./abs(tempX);
dy = tempY./abs(tempY);

% % update contributions
Dx = -dx; 
Dy = -dy;


        
