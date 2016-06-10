clear; close; clc;

finalTime   = 100; % simulation time
alpha       = 1.0;  % repulsion  distance 
rho         = 6.0;  % attraction distance 
w           = 0.5;  % weight factor (direction)
theta       = 2.0;  % angle threshold
s           = 0.5;  % speed constant
dt          = 0.1;  % time step
L           = 10.0; % boundary constraint (only if periodic)
g           = pi/2; % preferred direction 
pauseTime   = 0.1;  % pause time per animation
isAnime     = 0  ;  % animate results ? 1: ON, 0: OFF
isPeriodic  = 0  ;  % periodic boundaris ? 1: ON, 0: OFF


N_list     = [20; 50; 100; 200; 500];      % group size list
p_list     = (0:0.1:1.0)';                 % proportion list

% repetitions
numReps = 10; % number of repetitions

% initialize elongation
elong = zeros(numReps, length(p_list), length(N_list));

% parallel/serial version ?  (uncomment to use)
workingVersion();

tic % time start 
disp('simulation start...');
% start iterating
parfor N_idx=1:length(N_list)      % size
    N = N_list(N_idx);
    % monitor outer progress
    disp(['step: ', num2str(N_idx), ' out of ', num2str(length(N_list))]);
    % set all in/as function (eventually) for 
    % efficient parallelization ! 
    elong(:, :, N_idx) = ...
                     parallelFunction( p_list, numReps,     ...
                                       finalTime, N, alpha, ... 
                                       rho, w, s, dt, g,    ...
                                       L, theta, pauseTime, ...
                                       isAnime, isPeriodic );
    
end                            % size
disp('...simulation end');
toc % time lapsed 


% session is over
delete(gcp);


