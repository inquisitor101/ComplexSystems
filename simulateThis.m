function out = simulateThis(finalTime, N)

% position
Cx = zeros(N, finalTime); % horizontal position 
Cy = zeros(N, finalTime); % vertical   position 

% direction vector 
Vx = zeros(N, finalTime); % horizontal direction 
Vy = zeros(N, finalTime); % vertical   direction

% desired travel direction 
Dx = zeros(N, finalTime); % horizontal orientation 
Dy = zeros(N, finalTime); % vertical   orientation

% initialize position  (randomly)
Cx(:, 1) = rand(size(Cx(:, 1)));
Cy(:, 1) = rand(size(Cy(:, 1)));
% initialize direction (randomly)
Vx(:, 1) = rand(size(Vx(:, 1)));
Vy(:, 1) = rand(size(Vy(:, 1)));

% preliminary time step (initial: t=1)
for i=1:N   % individuals
    dx = 0; dy = 0;  % reset
    for j=1:i-1  % unroll loop 
        tempX = Cx(j, 1) - Cx(i, 1);
        tempY = Cy(j, 1) - Cy(j, 1); 
        dx = dx + tempX/abs(tempX); % normalized X
        dy = dy + tempY/abs(tempY); % normalized Y
    end
    for j=i+1:N  % unroll loop 
        tempX = Cx(j, 1) - Cx(i, 1);
        tempY = Cy(j, 1) - Cy(j, 1); 
        dx = dx + tempX/abs(tempX); % normalized X
        dy = dy + tempY/abs(tempY); % normalized Y
    end
    % update contributions
    Dx(i, 1) = dx + Vx(1, 1)/abs(Vx(1,1)); 
    Dy(i, 1) = dy + Vy(1, 1)/abs(Vy(1,1)); 
end         % individuals

for t=1:finalTime-1  % time 
       
    for i=1:N        % individuals
        dx = 0; dy = 0;  % reset
        for j=1:i-1  % unroll loop 
            tempX = Cx(j, t+1) - Cx(i, t+1);
            tempY = Cy(j, t+1) - Cy(j, t+1); 
            dx = dx + tempX/abs(tempX); % normalized X
            dy = dy + tempY/abs(tempY); % normalized Y
        end
        for j=i+1:N  % unroll loop 
            tempX = Cx(j, t+1) - Cx(i, t+1);
            tempY = Cy(j, t+1) - Cy(j, t+1); 
            dx = dx + tempX/abs(tempX); % normalized X
            dy = dy + tempY/abs(tempY); % normalized Y
        end
        % update contributions
        Dx(i, t+1) = dx + Vx(1, 1)/abs(Vx(1,t+1)); 
        Dy(i, t+1) = dy + Vy(1, 1)/abs(Vy(1,t+1));
        
    end              % individuals
    
end                  % time