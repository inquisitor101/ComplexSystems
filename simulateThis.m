% function out = simulateThis(finalTime, alpha, rho,...
%                             s, dt, p, w, g, N, ...
%                             isAnime, pauseTime)

% % % % % % % % % % % % % % % % % % % % % % % % %
finalTime   = 100; % simulation time
alpha       = 0.5; % repulsion  distance 
rho         = 2.0; % attraction distance 
w           = 0.5; % weight factor
s           = 0.5; % speed constant
dt          = 0.1; % time step
g           = [0; 0];
N           = 30;
p           = 0.1;
maxInformed = N*p;
L           = 10;
pauseTime   = 0.1;
isAnime     = 1;
% % % % % % % % % % % % % % % % % % % % % % % % %

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
Cx(:, 1) = L*rand(N, 1);
Cy(:, 1) = L*rand(N, 1);
% initialize direction (randomly)
Vx(:, 1) = rand(N, 1);
Vy(:, 1) = rand(N, 1);

% max number of informed individuals ~ p
maxInformed = N*p;
 
for t=1:finalTime-1  % time 
    
    for i=1:N   % individuals
        
        A = getDistance(i, L, Cx(:, t), Cy(:, t), N);
        
        t1 = nanmin(A);
        dist = sort(t1)';       % sorted distance
        idx  = dist <= alpha;   % repeled indices
        idx2 = dist <= rho;     % attract indices
        
        % repel
        if sum(idx) ~= 0
            [Dx(idx, t+1), Dy(idx, t+1)] = ...
                repelNeighbors  (idx, i, ...
                                 Cx(:, t), Cy(:, t) );
        elseif sum(idx2) ~= 0
        % attract     
            [Dx(idx2, t+1), Dy(idx2, t+1)] = ...
                attractNeighbors(idx2, i, ...
                                 Cx(:, t), Cy(:, t), ...
                                 Vx(:, t), Vy(:, t) );
        end
        
        % convert d to d^chapeau
        % make sure not to divide by zero !! 
        if Dx(i, t+1) ~= 0 
            Dx(i, t+1) = Dx(i, t+1)/abs(Dx(i, t+1));
        end
        if Dy(i, t+1) ~= 0
            Dy(i, t+1) = Dy(i, t+1)/abs(Dy(i, t+1));
        end
        
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
    
    % update desired direction
    temp = atan2( Dy(:, t+1), Dx(:, t+1) );
    % guassian distribution
    randomRotation = gaussianDist([-1, +1], 0.01, N);
    temp = temp + randomRotation; 
    
    Dx(:, t+1) = cos(temp); % horizontal location
    Dy(:, t+1) = sin(temp); % vertical   location
    
    % update direction
    Dangle = atan2(Dy(:, t+1), Dx(:, t+1)); % desired   angle
    Vangle = atan2(Vy(:, t), Vx(:, t));     % direction angle
    theta  = atan2(Cy(:, t), Cx(:, t));     % distance  angle

    idx = (Vangle-Dangle) < theta*dt; % check condition
    sgn = sign(Vangle-Dangle);        % orientation sense
    % condition near ( difference < theta * dt ) 
    Vx(idx, t+1) = Dx(idx, t+1);
    Vy(idx, t+1) = Dy(idx, t+1);
    idx2 = imcomplement(idx);   % implement else condition
    % condition far  ( difference > theta * dt )
    Vx(idx2, t+1) = Vx(idx2, t) + sgn(idx2).*theta(idx2).*dt;
    Vy(idx2, t+1) = Vy(idx2, t) + sgn(idx2).*theta(idx2).*dt;
    
    % update position
    Cx(:, t+1) = Cx(:, t) + Vx(:, t+1).*s*dt;
    Cy(:, t+1) = Cy(:, t) + Vy(:, t+1).*s*dt;
    
    % apply periodic boundaries (torus-like)
    Cx(:, t+1) = mod( Cx(:, t+1), L); % horizontal boundary
    Cy(:, t+1) = mod( Cy(:, t+1), L); % vertical   boundary
    
    % plot
    if isAnime
        for i=1:maxInformed
           idx_x = abs( Cx(i, t+1) - Cx(i, t) );
           idx_y = abs( Cy(i, t+1) - Cy(i, t) );
           if idx_x < 0.5*L && idx_y < 0.5*L
               plot([Cx(i, t), Cx(i, t+1)], [Cy(i, t), Cy(i, t+1)], 'r-','markersize',4)
               axis([0 L 0 L]);
               hold on
               plot(Cx(i, t+1), Cy(i, t+1), 'r.', 'markersize', 10)
               xlabel('X position')
               ylabel('Y position')
               title(['time: ',num2str(t), ...
               '   informed (red): ',num2str(maxInformed)]);
           end
        end
        for i=maxInformed+1:N
           idx_x = abs( Cx(i, t+1) - Cx(i, t) );
           idx_y = abs( Cy(i, t+1) - Cy(i, t) );
           if idx_x < 0.5*L && idx_y < 0.5*L
               plot([Cx(i, t), Cx(i, t+1)], ...
                    [Cy(i, t), Cy(i, t+1)], ...
                    'b-','markersize',4)
               axis([0 L 0 L]);
               hold on
               plot(Cx(i, t+1), Cy(i, t+1), 'b.', 'markersize', 10)
               xlabel('X position')
               ylabel('Y position')
               title(['time: ',num2str(t), ...
               '   informed (red): ',num2str(maxInformed)]);
           end
        end
    end
   % let's check things out 
   if isAnime
       getframe();
       pause(pauseTime); hold off
   end
end                  % time

if isAnime
    clf(); % clear figure 
end

