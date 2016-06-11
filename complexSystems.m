clear; clc;
%close;

finalTime   = 100; % simulation time
alpha       = 1.0;  % repulsion  distance 
rho         = 6.0;  % attraction distance 
w           = 0.5;  % weight factor (direction)
theta       = 2.0;  % angle threshold
s           = 1.0;  % speed constant
dt          = 0.2;  % time step
L           = 1.0; % boundary constraint (only if periodic)
g           = pi/3; % preferred direction 
pauseTime   = 0.05;  % pause time per animation
isAnime     = 0  ;  % animate results ? 1: ON, 0: OFF
isPeriodic  = 0  ;  % periodic boundaris ? 1: ON, 0: OFF


N_list     = [10; 30; 50; 100];      % group size list
p_list     = (0.1:0.1:1.0)';                 % proportion list

% repetitions
numReps = 10; % number of repetitions

% initialize elongation
elong = zeros(numReps, length(p_list), length(N_list));
% initialize group direction
vec   = zeros(numReps, length(p_list), length(N_list));
% initialise accuracy
acc   = zeros(numReps, length(p_list), length(N_list));

% parallel/serial version ?  (uncomment to use)
workingVersion();

tic % time start 
disp('simulation start...');
% start iterating
for N_idx=1:length(N_list)      % size
    N = N_list(N_idx);
    % monitor outer progress
    disp(['step: ', num2str(N_idx), ' out of ', num2str(length(N_list))]);
    % set all in/as function (eventually) for 
    % efficient parallelization ! 
    [elong(:, :, N_idx), vec(:, :, N_idx), acc(:, :, N_idx)] = ...
                     parallelFunction( p_list, numReps,     ...
                                       finalTime, N, alpha, ... 
                                       rho, w, s, dt, g,    ...
                                       L, theta, pauseTime, ...
                                       isAnime, isPeriodic );
    
end                            % size
disp('...simulation end');
toc % time lapsed 

%% Plot average accuracy
%Sv = mean(vec,1) - g;
figure(1)
clf
Sa = 1 - mean(acc,1)/sqrt(2);
hold on
plot(p_list,Sa(1,:,1))
plot(p_list,Sa(1,:,2))
plot(p_list,Sa(1,:,3))
plot(p_list,Sa(1,:,4))
hold off

%% Plot average elongation
figure(2)
clf
Se = mean(elong,1);
plot(p_list,Se(1,:,3))

%% session is over
delete(gcp);


