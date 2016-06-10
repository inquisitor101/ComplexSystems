
clear; clc; close 
finalTime   = 100; % simulation time
alpha       = 1.0; % repulsion  distance 
rho         = 6.0; % attraction distance 
w           = 0.5; % weight factor
s           = 1.0; % speed constant
dt          = 0.2; % time step
g           = pi/2;
N           = 10;
p           = 0.2;
maxInformed = N*p;
L           = 10.0;
theta       = 2.0;
pauseTime   = 0.0;
isAnime     = 1;
isPeriodic  = 1;

% position
Cx = zeros(N, finalTime); Cy = zeros(N, finalTime); 
% direction vector 
Vx = zeros(N, finalTime); Vy = zeros(N, finalTime); 
% desired direction
Dx = zeros(N, finalTime); Dy = zeros(N, finalTime);
% centroid
Xc = zeros(finalTime,1);    Yc = zeros(finalTime,1);
% group direction 
h  = zeros(finalTime,1);

% initialize position  (randomly)
Cx(:, 1) = 0.25*L*rand(N, 1)+0.375*L; % centered in L-by-L
Cy(:, 1) = 0.25*L*rand(N, 1)+0.375*L; % centered in L-by-L
% initialize direction (randomly)
Vx(:, 1) = 2*pi*rand(N, 1);
Vy(:, 1) = 2*pi*rand(N, 1);
% initialize centroids
Xc(1) = mean(Cx(:, 1)); Yc(1) = mean(Cy(:, 1));

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% step 1: traverse time
for t=1:finalTime-1
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 2: sweep through particles
    for i=1:maxInformed
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.1: distance
        if isPeriodic
            A = getDistance(i, L, Cx(:, t), Cy(:, t), N);
            dist = nanmin(A);
        else
            xR = Cx(i, t) - Cx(:, t); yR = Cy(i, t) - Cy(:, t);
            xR(i) = 2*rho; yR(i) = 2*rho; % exclude self (no delete)
            dist = sqrt(xR.^2 + CyR.^2);
        end
        idx  = dist <= alpha;   % repeled indices
        idx2 = dist <= rho;     % attract indices
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.2: neighbors
        if sum(idx) ~= 0       % repel
            
            [Dx(i, t+1), Dy(i, t+1)] = repelNeighbors  (idx, i, Cx(:, t), Cy(:, t) );
        
        elseif sum(idx2) ~= 0  % attract
            
            [Dx(i, t+1), Dy(i, t+1)] = attractNeighbors(idx2, i, Cx(:, t), Cy(:, t),  Vx(:, t), Vy(:, t) );
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.3: normalize 
        D = atan2(Dy(i, t+1), Dx(i, t+1) ); % convert to angle
        D = D/(2*pi);  % normalize using 2*pi
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.4: normalize using preferred direction
        D = D + w*g;
        D = D/(2*pi);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.5: update Dx Dy
        Dx(i, t+1) = cos(D); Dy(i, t+1) = sin(D);
    end
    for i=maxInformed+1:N
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.6: distance
        if isPeriodic
            A = getDistance(i, L, Cx(:, t), Cy(:, t), N);
            dist = nanmin(A);
        else
            xR = Cx(i, t) - Cx(:, t); yR = Cy(i, t) - Cy(:, t);
            xR(i) = 2*rho; yR(i) = 2*rho; % exclude self (no delete)
            dist = sqrt(xR.^2 + CyR.^2);
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.7: neighbors
        if sum(idx) ~= 0       % repel
            
            [Dx(i, t+1), Dy(i, t+1)] = repelNeighbors  (idx, i, Cx(:, t), Cy(:, t) );
        
        elseif sum(idx2) ~= 0  % attract
            
            [Dx(i, t+1), Dy(i, t+1)] = attractNeighbors(idx2, i, Cx(:, t), Cy(:, t),  Vx(:, t), Vy(:, t) );
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.8: normalize 
        D = atan2(Dy(i, t+1), Dx(i, t+1) ); % convert to angle
        D = D/(2*pi);  % normalize using 2*pi 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % step 2.9: update Dx Dy
        Dx(i, t+1) = cos(D); Dy(i, t+1) = sin(D);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 3: add noise + update Dx Dy
    angle = atan2( Dy(:, t+1), Dx(:, t+1) );
    noise = pi/2*gaussianDist([-1, +1], 0.01, N);
    angle = angle + noise;
    Dx(:, t+1) = cos(angle); Dy(:, t+1) = sin(angle);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 4: update Vx Vy
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 4.1: get Vangle Dangle
    Dangle = atan2(Dy(:, t+1), Dx(:, t+1)); % desired   angle
    Vangle = atan2(Vy(:, t),   Vx(:, t));   % direction angle
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 4.2: get indices satisfying condition: theta*dt
    idx = abs(Dangle-Vangle) < theta*dt; 
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 4.3: update according to criteria
    Vx(idx, t+1) = Dx(idx, t+1); Vy(idx, t+1) = Dy(idx, t+1);
    idx2 = imcomplement(idx); sgn = sign(Dangle-Vangle);
    Vangle(idx2) = Vangle(idx2) + sgn(idx2)*theta*dt; 
    Vx(idx2, t+1) = cos(Vangle(idx2)); Vy(idx2, t+1) = sin(Vangle(idx2));
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 5: update position Cx Cy
    Cx(:, t+1) = Cx(:, t) + Vx(:, t+1)*s*dt;
    Cy(:, t+1) = Cy(:, t) + Vy(:, t+1)*s*dt;
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 6: get centroid Xc Yc
    Xc(t+1) = mean(Cx(:, t+1)); Yc(t+1) = mean(Cy(:, t+1));
    
    % this part is under testing...
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 7: group direction
    h(t+1) = atan2(Yc(t+1)-Yc(t), Xc(t+1)-Xc(t) );
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 8: apply periodic boundaries 
    if isPeriodic
        Cx(:, t+1) = mod( Cx(:, t+1), L); % horizontal boundary
        Cy(:, t+1) = mod( Cy(:, t+1), L); % vertical   boundary
    end
    
    
    if isAnime
       animateThis(maxInformed, N, L, t, h, Cx, Cy, Xc, Yc, pauseTime); 
    end
end

clf()