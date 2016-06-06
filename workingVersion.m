function [] = workingVersion()

prompt     = '\nPlease enter which mode to use: \n<0>: serial\n<1>: parallel\n\n';
isParallel = input(sprintf(prompt));
clc;
Nproc      = gcp('nocreate');
if isParallel == 1 
    % parallel initialization
    if isempty(Nproc)  % parallel
        parpool        % begin parallel session
    end
    disp('Parallel Version ON');
else 
    % serial initialization
    if ~isempty(Nproc) % check parallel
        delete(gcp);   % delete current session
    end
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
    disp('Serial Verion ON');
end