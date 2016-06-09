function elong = boundingBox(Cx, Cy, h)

% direction of travel (w.r.t. centroid)
groupDirection    = h;
% direction of travel of each particle/individual
particleDirection = atan2(Cy, Cx);

% Tolerance ? 
TOL = 1e-2;
% width of elongation (parallel length)
par = abs(particleDirection - groupDirection);
par1 = par(find(par<=TOL & par >= -TOL) );

par2 = par(find(par<= (pi+TOL) & par>= (par-TOL)) );
par = par1 + par2;

par_min = min(par); par_max = max(par);

% hight of elongation (perpendicular   length)
perp = abs(particleDirection - groupDirection);
perp = perp(find(perp <= (pi/2+TOL) & perp >= (-pi/2-TOL)) );

perp_min = min(par); perp_max = max(par);
% elongation (parallel length over perpendicular length)
elong = (par_max - par_max) / (perp_max - perp_min);



