function out = simulateThis(finalTime, alpha, rho, N)

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
Cx(:, 1) = rand(size( Cx(:, 1)) );
Cy(:, 1) = rand(size( Cy(:, 1)) );
% initialize direction (randomly)
Vx(:, 1) = rand(size( Vx(:, 1)) );
Vy(:, 1) = rand(size( Vy(:, 1)) );



for t=1:finalTime-1  % time 
    
    
    for i=1:N
        % distance check
        dist = getDistance(i, Cx, Cy, N);

        % repel
        if sum(dist < alpha) ~= 0   % within alpha range

            [Dx(i, t+1), Dy(i, t+1)] = ...
                repelNeighbors  (N, i, ...
                                 Cx(:, t+1), Cy(:, t+1) );
        % attract
        elseif sum(dist < rho) ~= 0 % within rho interaction 
            
            [Dx(i, t+1), Dy(i, t+1)] = ...
                attractNeighbors(N, i, ...
                                 Cx(:, t+1), Cy(:, t+1), ...
                                 Vx(:, t+1), Vy(:, t+1) );

        end
    end
    
end                  % time

