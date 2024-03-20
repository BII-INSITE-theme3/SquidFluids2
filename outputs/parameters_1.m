% Parameters needed to run the 2D flow code

% Parameters:

% Stokeslet parameters
rho = 20; % Number of stokeslets per unit length.
eps = 0.01; % Regularisaton parameter.

% Channel geometry
Lt = 10; % Length of top segment.
Lm = 8; % Length of transition region.
Lb = 10; % Length of bottom segment.
theta = pi/4; % Angle of right transition region to horizontal.
Ltot = Lt+sin(pi/4)*Lm+Lb; % Total height of system simulated.
Ptx = 10; % Position of top point of right boundary.
Pty = Lt+Lm/2; % Position of top point of right boundary.

% Appendage geometry
dsep = 1; % Appendage separation.
psi = 0; % Angle of inclination between appendage pairs (Rad).
PRAx = 5; % Position of right appendage in x.
PRAy = 5; % Position of right appendage in y.

% Flow parameters
U0 = -1; % Background flow strength Max.

% Underlying space parameters.
nptx = 100; % Solver points in x direction.
npty = nptx; % Solver points in y direction.
x = linspace(-(Ptx+1)-5,(Ptx+1)+5,nptx); % Solver x coords.
y = linspace(-(Ptx+1)-5,(Ptx+1)+5,npty); % Solver y coords.
