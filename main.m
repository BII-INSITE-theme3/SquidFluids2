% Title: Code to get the fluid flow from a set of Regularized Stokeslets.
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

%% Output for user and begin the execution timer
disp(['Main running correctly.']);
tic;

%% Add the function files need to run
addpath('functions/') 

%% Set parameters
parameters % Set the parameters

%% Determine if the parameters have been pre-run
[pre_found,file_index,stks,F,dUflowy,Uflowx,Uflowy] = check_existing_runs;

%% Calculation loop
% If there is no pre-existing run, calculate the values
if ~pre_found

    disp('No pre-existing run data exists. Now calculating the flow profile.')

    %% Initialise parallelisation pool
    parpool; % Intialise the parallel workers % (Optional, to speed up runtime)
    
    %% Set channel geometry
    % Set the system geometry
    stks = geometry(rho,Lt,Lm,Lb,theta,Ptx,Pty,dsep,psi,PRAx,PRAy);
    
    disp('Geometry initialised.')
    
    %% Solve for the forces
    F = getForces(stks,eps,U0);
    
    disp('Stokeslet forces caclulated.')
    
    %% Optional: Check the boundary flow is as expected
    % Calculates the vales of the resulting flow from the forces at the points of the Stokeslets. 
    % This allows us to show that the re-constituted flow field is exactly the boundary conditions
    
    a = find(stks(:,3)==1);
    [~,Uflowy] = calculateFlowPath(stks,F,stks(a,1),stks(a,2),eps);
    dUflowy = mean(Uflowy);
    % n=1; % Plot density (to make visualising easier for large numbers of points)
    % quiver(stks(:,1),stks(:,2),Uflowx(1:n:end),Uflowy(1:n:end)) % Plot the vector field
    
    disp('Correction flow caclulated.')
    
    %% Get the flow over the space
    [Uflowx,Uflowy] = calculateFlowGrid(stks,F,x,y);
    Uflowy = Uflowy - dUflowy;

    disp('Full system flow grid caclulated.')

    %% Optional: Save the outputs
    % Check the output folder exists
    if ~exist('outputs','file') == 7
        mkdir outputs
    end
    
    %% Save the outputs from the run
    save(['outputs/output_' num2str(file_index)], 'stks','F','dUflowy','Uflowx','Uflowy'); % Save the calculated values.
    sourceFile = 'parameters.m'; % Current parameter file
    destinationFile = ['outputs/parameters_' num2str(file_index) '.m']; % Copy path
    
    % Copy the parameter file
    copyfile(sourceFile, destinationFile);

end

%% Plot output flow field
% This section calculates the local magnitude of the flow the create the colourised background. 
% Then the quiver plot can be overlaid on top of this to give the directions of the flow.

n=15; % Plot coarseness
Umag = sqrt(Uflowx.^2 + Uflowy.^2); % Get flow field magnitude
imagesc(x,y,Umag) % Plot local flow strenght as a background
hold on
scatter(stks(:,2),stks(:,1),2,'r') % Plot the Stokeslets
quiver(x(1:n:end),y(1:n:end),Uflowy(1:n:end,1:n:end),Uflowx(1:n:end,1:n:end),2,'Color','y') % Plot the vector field
axis equal

%% Close the pool of parallel workers
delete(gcp('nocreate')) % (Optional to speed up runtime)

%% Finalize the run
time = toc;
disp(['Code execution time: ' num2str(time) 's'])
