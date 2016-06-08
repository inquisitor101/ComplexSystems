function [Dx, Dy] = attractNeighbors(idx2, i, Cx, Cy, Vx, Vy)


tempX = Cx(idx2) - Cx(i);
tempY = Cy(idx2) - Cy(i);
dx = tempX./abs(tempX);
dy = tempY./abs(tempY);
% correct for case division by zero 
dx(isnan(dx)) = 0;
dy(isnan(dy)) = 0;

tempX = Vx./abs(Vx);
tempY = Vy./abs(Vy);

dv_x = sum(tempX);
dv_y = sum(tempY);

Dx = dx + dv_x;
Dy = dy + dv_y;



