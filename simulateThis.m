function out = simulateThis(finalTime, alpha, rho,...
                            p, w, g, N, maxInformed)

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
    
    for i=1:N   % individuals
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

        end  % repel/attract condition
        
        % convert d to d^chapeau
        Dx(i, t+1) = Dx(i, t+1)/abs(Dx(i, t+1));
        Dy(i, t+1) = Dy(i, t+1)/abs(Dy(i, t+1)); 
        
        % update informed individuals
        if i <= maxInformed 
            % assume first 'maxInformed' are the 
            % informed individuals for simplicity 
            % and book-keeping purposes.
            tempX = Dx(i, t+1) + w*g(1);
            tempX = tempX/abs(tempX);
            tempY = Dy(i, t+1) + w*g(2);
            tempY = tempY/abs(tempY);
            Dx(i, t+1) = tempX;
            Dy(i, t+1) = tempY;
        end
        
    end         % individuals
    
    % update direction
    temp = atan2( Dy(:, t+1), Dx(:, t+1) );
    %
    % TODO: 
    % figure how to get a gaussian distribution
    % about 0, standDeviation = 0.01
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    temp = temp + randomRotation; % guassian distribution
    Dx(:, t+1) = cos(temp);
    Dy(:, t+1) = sin(temp);
    
end                  % time






