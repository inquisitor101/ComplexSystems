clear; close; clc;

finalTime   = 10 ; % simulation time
alpha       = 0.5; % repulsion  distance 
rho         = 2.0; % attraction distance 
w           = 0.5; % weight factor
s           = 0.5; % speed constant
dt          = 0.1; % time step
maxInformed = 10 ; % number of informed individuals

% g(1): horizontal, g(2): vertical
g          = [0; 1];       % preferred direction 

N_list     = (1:10)';      % group size list
p_list     = (0:0.1:1.0)'; % proportion list
% repetitions
numReps = 10; % number of repetitions

% parallel/serial version ?  (uncomment to use)
% workingVersion();


tic % time start 
disp('simulation start...');
% start iterating
parfor N=1:length(N_list)      % size
    
    % set all in function eventually for 
    % efficient parallelization ! 
    for p_idx=1:length(p_list) % proportion
        p = p_list(p_idx);
        for r=1:numReps        % repetitions     
            simulateThis(finalTime, alpha, rho,...
                         p, w, g, N, maxInformed);
        end                    % repetitions
    end                        % proportion
    
end                            % size
disp('...simulation end');
toc % time lapsed 



% session is over
delete(gcp);


