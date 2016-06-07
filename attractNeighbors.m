function [Dx, Dy] = attractNeighbors(idx2, i, Cx, Cy, Vx, Vy)


tempX = Cx(idx2) - Cx(i);
tempY = Cy(idx2) - Cy(i);
dx = tempX./abs(tempX);
dy = tempY./abs(tempY);
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% TODO: 
% Vx(1)/abs(Vx(1))    ... doesn't make sense 
% double check if its always index = 1 1?!?!?!?
Dx = dx + Vx(1)/abs(Vx(1));
Dy = dy + Vy(1)/abs(Vy(1));



