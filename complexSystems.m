clear; close; clc;

groupSize  = 10;  % group size
finalTime  = 10;  % simulation time
w          = 0.5; % weight factor
s          = 0.5; % speed constant
dt         = 0.1; % time step

% position
Cx = zeros(groupSize, finalTime); % horizontal position x-axis
Cy = zeros(groupSize, finalTime); % vertical   position y-axis

% distance 
Dx = zeros(groupSize, finalTime); % horizontal distance 
Dy = zeros(groupSize, finalTime); % vertical   distance

% repetitions
numReps = 50; % # of repetitions

% parallel/serial version ?  (uncomment to use)
% workingVersion();

% start iterating
parfor N=1:groupSize  % individuals
    
    
    
end  % individuals



% session is over
delete(gcp);