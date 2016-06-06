clear; close; clc;

N_list     = 10;  % group size
finalTime  = 10;  % simulation time
w          = 0.5; % weight factor
s          = 0.5; % speed constant
dt         = 0.1; % time step

% repetitions
numReps = 50; % number of repetitions

% parallel/serial version ?  (uncomment to use)
% workingVersion();


tic % time start 
disp('simulation start...');
% start iterating
parfor N=1:N_list          % size
    % set all in function eventually for parallelization ! 
    
    for r=1:numReps     % repetitions     
        simulateThis(finalTime, N);
    end                 % repetitions
    
end                        % size
disp('...simulation end');
toc % time lapsed 

% session is over
delete(gcp);




