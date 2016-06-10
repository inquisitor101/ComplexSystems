function [box, elong] = boundingBox(Cx, Cy, theta)

% 2-vector matrix (to-be rotated)
matrixIN = [Cx'; Cy'];
% bounding box initialization
box      = zeros(2, 4);

% rotation matrix operator
Roperator = @(theta)[ cos(theta), -sin(theta); 
                      sin(theta),  cos(theta) ];
          
% reverse rotation (step: 1/2)  
R = Roperator(-theta);

% rotate according to rotation matrix
matrixOUT = (R*matrixIN)';

% vector decomposition:
%                       x | y
%                      ---|---
%                       # | #
% 
% % % % % % % % % % % % % % % % % % % %
M_min = min(matrixOUT)';
M_max = max(matrixOUT)';


% complement rotation (step: 2/2)
R = R';
  
box(:,1) = R*[M_min(1); M_min(2)];
box(:,2) = R*[M_min(1); M_max(2)];
box(:,3) = R*[M_max(1); M_max(2)];
box(:,4) = R*[M_max(1); M_min(2)];


elong = (M_max(1) - M_min(1)) /...
        (M_max(2) - M_min(2)) ;











