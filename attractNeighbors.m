function [Dx, Dy] = attractNeighbors(N, i, Cx, Cy, Vx, Vy)


dx = 0; dy = 0;  % reset
for j=1:i-1  % unroll loop 
    tempX = Cx(j) - Cx(i);
    tempY = Cy(j) - Cy(j); 
    dx = dx + tempX/abs(tempX); % normalized X
    dy = dy + tempY/abs(tempY); % normalized Y
end
for j=i+1:N  % unroll loop 
    tempX = Cx(j) - Cx(i);
    tempY = Cy(j) - Cy(j); 
    dx = dx + tempX/abs(tempX); % normalized X
    dy = dy + tempY/abs(tempY); % normalized Y
end

% update contributions
Dx = dx + Vx(1, 1)/abs(Vx(1)); 
Dy = dy + Vy(1, 1)/abs(Vy(1));

