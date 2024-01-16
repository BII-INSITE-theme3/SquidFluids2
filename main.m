% Title: Second attempt to do the inverse problem.
% Author: Stephen Williams.
% Notes: 1. Code is non-dimensionalised. 
% 1.1 The characteristic length is equal to the appendage radius. 
% 1.2 The characteristic time is the time it would take for the flow at the maximum on the entry point to move 1 appendage radius.
% 1.3 The viscocity is set so that the prefactor in the stokeslet is 1.
% 2. The boundaries are now arc-length parameterisation (instead of uniform
% along y). This improves the stokeslet distribution on the non-vertical
% portions of the boundaries.
%--------------------------------------------%

close all
clear all

%% Set parameters
parameters % Set the parameters, TBD

%% Initialise parallelisation pool
%parpool % Intialise the parallel workers % (Optional, to speed up runtime)

%% Set channel geometry
% Set the system geometry
stks = geometry(rho,Lt,Lm,Lb,theta,Ptx,Pty,dsep,psi,PRAx,PRAy);

%% Solve for the forces
F = getForces(stks,eps,U0);

%% Check the boundary flow is as expected
% Calculates the vales of the resulting flow from the forces at the points of the Stokeslets. This allows us to show that the re-constituted flow field is exactly the boundary conditions

% n=1; % Plot density (to make visualising easier for large numbers of points)
% [Uflowx,Uflowy] = calculateFlowPath(stks,F,stks(:,1),stks(:,2),eps);
% quiver(stks(:,1),stks(:,2),Uflowx(1:n:end),Uflowy(1:n:end)) % Plot the vector field

%% Get the flow over the space
 [Uflowx,Uflowy] = calculateFlowGrid(stks,F,x,y);

%% Plot output flow field
% This section calculates the local magnitude of the flow the create the colourised background. 
% Then the quiver plot can be overlaid on top of this to give the directions of the flow.

n=1; % Plot coarseness
Umag = sqrt(Uflowx.^2 + Uflowy.^2); % Get flow field magnitude
imagesc(x,y,Umag') % Plot local flow strenght as a background
hold on
scatter(stks(:,1),stks(:,2),2,'r') % Plot the Stokeslets
%quiver(x(1:n:end),y(1:n:end),Uflowx(1:n:end,1:n:end),Uflowy(1:n:end,1:n:end)) % Plot the vector field
axis equal

%% Close the pool of parallel workers
% delete(gcp('nocreate')) % (Optional to speed up runtime)
