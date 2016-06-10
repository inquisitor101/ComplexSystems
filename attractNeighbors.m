function [Dx, Dy] = attractNeighbors(idx2, i, Cx, Cy, Vx, Vy)

idx2(i) = 0; % make sure exclude self 

tempX = Cx(idx2) - Cx(i);
tempY = Cy(idx2) - Cy(i);
dx = tempX./sqrt(tempX.^2 + tempY.^2);
dy = tempY./sqrt(tempX.^2 + tempY.^2);
% correct for case division by zero 
dx(isnan(dx)) = 0;
dy(isnan(dy)) = 0;

dv_x = Vx./sqrt(Vx.^2 + Vy.^2);
dv_y = Vy./sqrt(Vx.^2 + Vy.^2);

Dx = sum(dx) + sum(dv_x);
Dy = sum(dy) + sum(dv_y);


