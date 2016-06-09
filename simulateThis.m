% function out = simulateThis(finalTime, alpha, rho,...
%                             s, dt, p, w, g, N, ...
%                             isAnime, pauseTime)

% % % % % % % % % % % % % % % % % % % % % % % % %
clear; close; clc;
finalTime   = 100; % simulation time
alpha       = 0.5; % repulsion  distance 
rho         = 2.0; % attraction distance 
w           = 1.0; % weight factor
s           = 1.0; % speed constant
dt          = 0.1; % time step
g           = [1; 0];
N           = 10;
p           = 0.1;
maxInformed = N*p;
L           = 2;
pauseTime   = 0.05;
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

% centroid location of the whole group (N)
groupCentroidX = zeros(finalTime);
groupCentroidY = zeros(finalTime);

% initialize position  (randomly)
Cx(:, 1) = L*rand(N, 1);
Cy(:, 1) = L*rand(N, 1);
% initialize direction (randomly)
Vx(:, 1) = rand(N, 1);
Vy(:, 1) = rand(N, 1);

% initial step (t = 1)
groupCentroidX = mean(Cx(:, 1));
groupCentroidY = mean(Cy(:, 1));

% max number of informed individuals ~ p
maxInformed = N*p;
 
% check p condition
if p > 1 || p < 0
    error('proportion must obey condition: 0.0 <= p <= 1.0');
end

for t=1:finalTime-1  % time 
    
    for i=1:N   % individuals
        
        A = getDistance(i, L, Cx(:, t), Cy(:, t), N);
        
        dist = nanmin(A);
        idx  = dist <= alpha;   % repeled indices
        idx2 = dist <= rho;     % attract indices
        
        % repel
        if sum(idx) ~= 0
            [Dx(i, t+1), Dy(i, t+1)] = ...
                repelNeighbors  (idx, i, ...
                                 Cx(:, t), Cy(:, t) );
        elseif sum(idx2) ~= 0
        % attract     
            [Dx(i, t+1), Dy(i, t+1)] = ...
                attractNeighbors(idx2, i, ...
                                 Cx(:, t), Cy(:, t), ...
                                 Vx(:, t), Vy(:, t) );
        end
        
        % convert d to d^chapeau
        % make sure not to divide by zero !!
        d = sqrt(Dx(i, t+1).^2 + Dy(i, t+1).^2);
        if d ~= 0 
            Dx(i, t+1) = Dx(i, t+1)/d;
            Dy(i, t+1) = Dy(i, t+1)/d;
        end
        
        % update informed individuals
        if i <= maxInformed 
            % assume first 'maxInformed' are the 
            % informed individuals for simplicity 
            % and book-keeping purposes.
            
            tempX = Dx(i, t+1) + w*g(1);
            tempY = Dy(i, t+1) + w*g(2);
            d = sqrt(tempX.^2 + tempY.^2);
            % zero division condition
            if d ~= 0
                Dx(i, t+1) = tempX/d;
                Dy(i, t+1) = tempY/d;
            end
        end
        
    end         % individuals
    
    % update desired direction
    temp = atan2( Dy(:, t+1), Dx(:, t+1) );
    % guassian distribution
    randomRotation = gaussianDist([-1, +1], 0.01, N);
    temp = temp + randomRotation; 
    
    Dx(:, t+1) = cos(temp); % horizontal location
    Dy(:, t+1) = sin(temp); % vertical   location
    
    % update direction: Vx and Vy
    % step 1: get angles 
    Dangle = atan2(Dy(:, t+1), Dx(:, t+1)); % desired   angle
    Vangle = atan2(Vy(:, t), Vx(:, t));     % direction angle
    theta  = atan2(Cy(:, t), Cx(:, t));     % distance  angle
    % step 2: get indices of difference criteria
    idx = abs(Vangle-Dangle) < theta*dt; % check condition
    sgn = sign(Vangle-Dangle);           % orientation sense
    % step 3: update accordingly !
    % condition near ( difference < theta * dt ) 
    Vx(idx, t+1) = Dx(idx, t+1);
    Vy(idx, t+1) = Dy(idx, t+1);
    idx2 = imcomplement(idx);   % implement else condition
    % condition far  ( difference > theta * dt )
    Vx(idx2, t+1) = Vx(idx2, t) + sgn(idx2).*theta(idx2)*dt;
    Vy(idx2, t+1) = Vy(idx2, t) + sgn(idx2).*theta(idx2)*dt;
    
    % update position
    Cx(:, t+1) = Cx(:, t) + Vx(:, t+1)*s*dt;
    Cy(:, t+1) = Cy(:, t) + Vy(:, t+1)*s*dt;
    
    % apply periodic boundaries (torus-like)
    Cx(:, t+1) = mod( Cx(:, t+1), L); % horizontal boundary
    Cy(:, t+1) = mod( Cy(:, t+1), L); % vertical   boundary
    
    % centroid 
    groupCentroidX(t+1) = mean(Cx(:, t+1)); % horizontal center
    groupCentroidY(t+1) = mean(Cy(:, t+1)); % vertical   center
    
    % plot
    if isAnime
        for i=1:maxInformed
           idx_x = abs( Cx(i, t+1) - Cx(i, t) );
           idx_y = abs( Cy(i, t+1) - Cy(i, t) );
           if idx_x < 0.5*L && idx_y < 0.5*L
               plot([Cx(i, t), Cx(i, t+1)], ...
                    [Cy(i, t), Cy(i, t+1)], ...
                    'r-','markersize',4)
               axis([0 L 0 L]);
               hold on
               plot(Cx(i, t+1), Cy(i, t+1), 'r.', ...
                       'markersize', 10)
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
               plot(Cx(i, t+1), Cy(i, t+1), 'b.', ...
                    'markersize', 10)
               xlabel('X position')
               ylabel('Y position')
               title(['time: ',num2str(t), ...
               '   informed (red): ',num2str(maxInformed)]);
           end
        end
    end
   % let's check things out 
   if isAnime
       % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
       % debugging here (remove) ------------------------------
       weird = sum(isnan(Cx(:, t+1)));
       legend(['# NaN: ',num2str(weird)]);    
       % debugging here (remove) ------------------------------
       % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
       getframe();
       pause(pauseTime); hold off
   end
end                  % time

if isAnime
    clf(); % clear figure 
end

