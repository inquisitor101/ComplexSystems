clear; close; clc;

finalTime   = 100; % simulation time
alpha       = 0.5; % repulsion  distance 
rho         = 2.0; % attraction distance 
w           = 0.5; % weight factor
s           = 0.5; % speed constant
dt          = 0.1; % time step
pauseTime   = 0.1; % pause time per animation
isAnime     = 0  ; % animate results ? 1: ON, 0: OFF

% direction preference, properties:
% g(1): horizontal  (x-axis) 
% g(2): vertical    (y-axis)
%  (no change)  0 <= g(#) <= 1  (change direction) 
g          = [0; 1];       % preferred direction 

N_list     = (1:10)';      % group size list
p_list     = (0:0.1:1.0)'; % proportion list

% repetitions
numReps = 1; % number of repetitions

% parallel/serial version ?  (uncomment to use)
% workingVersion();


tic % time start 
disp('simulation start...');
% start iterating
parfor N=1:length(N_list)      % size
    
    % set all in/as function (eventually) for 
    % efficient parallelization ! 
    for p_idx=1:length(p_list) % proportion
        p = p_list(p_idx);
        for r=1:numReps        % repetitions     
            simulateThis(finalTime, alpha, rho,...
                         s, dt, p, w, g, N, ...
                         isAnime, pauseTime);
        end                    % repetitions
    end                        % proportion
    
end                            % size
disp('...simulation end');
toc % time lapsed 


% session is over
delete(gcp);


