function [elong, vec, acc] = simulation(finalTime, N, alpha, rho, w, s, dt, ...
                              g, p, L, theta, pauseTime, isAnime, ...
                              isPeriodic)


% total informed population
maxInformed = round(N*p);
% position
Cx = zeros(N, finalTime); Cy = zeros(N, finalTime); 
% direction vector 
Vx = zeros(N, finalTime); Vy = zeros(N, finalTime); 
% desired direction
Dx = zeros(N, finalTime); Dy = zeros(N, finalTime);
% centroid
Xc = zeros(finalTime,1);  Yc = zeros(finalTime,1);
% group direction 
h  = zeros(finalTime,1);

% initialize position  (randomly)
Cx(:, 1) = 0.25*L*rand(N, 1)+0.375*L; % centered in L-by-L
Cy(:, 1) = 0.25*L*rand(N, 1)+0.375*L; % centered in L-by-L
% initialize direction (randomly)
Vx(:, 1) = rand(N, 1);
Vy(:, 1) = rand(N, 1);
% initialize centroids
Xc(1) = mean(Cx(:, 1)); Yc(1) = mean(Cy(:, 1));

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% step 1: traverse time
for t=1:finalTime-1
%     % monitor progress
%     disp(['step: ', num2str(t+1), ' out of ', num2str(finalTime)]);
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
            dist = sqrt(xR.^2 + yR.^2);
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
        Dx(i, t+1) = cos(D); Dy(i, t+1) = sin(D);
        %D = D/sqrt(Dx(i, t+1)^2 + Dy(i, t+1)^2);  % normalize using 2*pi
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.4: normalize using preferred direction
        Dx = Dx + w*cos(g);
        Dy = Dy + w*sin(g);
        D = atan2(Dy(i, t+1), Dx(i, t+1) ); % convert to angle
        Dx(i, t+1) = cos(D); Dy(i, t+1) = sin(D);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.5: update Dx Dy
        %Dx(i, t+1) = cos(D); Dy(i, t+1) = sin(D);
    end
    for i=maxInformed+1:N
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % step 2.6: distance 
        if isPeriodic
            A = getDistance(i, L, Cx(:, t), Cy(:, t), N);
            dist = nanmin(A);
        else
            xR = Cx(:, t) - Cx(i, t); yR = Cy(:, t) - Cy(i, t);
            xR(i) = 2*rho; yR(i) = 2*rho; % exclude self (no delete)
            dist = sqrt(xR.^2 + yR.^2);
        end
        idx  = dist <= alpha;   % repeled indices
        idx2 = dist <= rho;     % attract indices
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
%         D = D/sqrt(Dx(i, t+1)^2 + Dy(i, t+1)^2);  % normalize using 2*pi 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % step 2.9: update Dx Dy
        Dx(i, t+1) = cos(D); Dy(i, t+1) = sin(D);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 3: add noise + update Dx Dy
%     angle = atan2( Dy(:, t+1), Dx(:, t+1) );
%     noise = pi/2*gaussianDist([-1, +1], 0.01, N);
%     angle = angle + noise;
%     Dx(:, t+1) = cos(angle); Dy(:, t+1) = sin(angle);
    
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
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % step 9: elongation
    [box, e] = boundingBox(Cx(:, t+1), Cy(:, t+1), h(t+1));
    
    if isAnime
       animateThis(maxInformed, N, L, t, h, Cx, Cy, Xc, Yc, pauseTime, isPeriodic, box, e); 
    end
end

if isAnime
    clf()
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% step 10: elongation
[~, elong] = boundingBox(Cx(:, end), Cy(:, end), h(end));

% condition check
if finalTime < 51
    error('max time must be at leasat 51 time steps. Change !');
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% step 11: group direction as defined in paper
rise = Yc(finalTime) - Yc(finalTime-50);
run  = Xc(finalTime) - Xc(finalTime-50);

vec  = atan2(rise,run);

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% step 12: group accuracy
acc = sqrt(sum(2*(1-cos(atan2(Vy(:,end),Vx(:,end)) - g)))/(2*N));

