function elong = boundingBox(Cx, Cy, h)

% definitions
par_min  = zeros(2, 1); par_max  = zeros(2, 1);
perp_min = zeros(2, 1); perp_max = zeros(2, 1);
% direction of travel (w.r.t. centroid)
groupDirection    = h;
% direction of travel of each particle/individual
particleDirection = atan2(Cy, Cx);

% Tolerance, angle within +/- 14.4 degrees of direction 
TOL = 0.04*(2*pi);

% angular difference
dif = abs(particleDirection - groupDirection);
% width of elongation (parallel length)
par1 = dif < (0.0 + TOL) & dif > (0.0 - TOL);
par2 = dif < (pi  + TOL) & dif > (pi  - TOL);
par  = par1 | par2;

par_min(1) = min(Cx(par)); par_max(1) = max(Cx(par)); % x
par_min(2) = min(Cy(par)); par_max(2) = max(Cy(par)); % y

% hight of elongation (perpendicular length)
perp1 = dif < (pi/2   + TOL) & dif > (pi/2   - TOL);
perp2 = dif < (3*pi/2 + TOL) & dif > (3*pi/2 - TOL);
perp  = perp1 | perp2;

if sum(perp) ~= 0
    perp_min(1) = min(Cx(perp)); perp_max(1) = max(Cx(perp)); % x
    perp_min(2) = min(Cy(perp)); perp_max(2) = max(Cy(perp)); % y
else
    perp_min(:) = 0; perp_max(:) = 0;
end
% distance calculation
xR = perp_max(1) - perp_min(1);
yR = perp_max(2) - perp_min(2);
rise = sqrt( xR*xR + yR*yR ); % case: perpendicular 

xR = par_max(1) - par_min(1);
yR = par_max(2) - par_min(2);
run  = sqrt( xR*xR + yR*yR ); % case: parallel/collinear

% elongation (parallel length over perpendicular length)
elong = rise/run;



