% Title: Streamline calculations.
% Author: Stephen Williams.
%--------------------------------------------%

clear all
close all

%% Load in the data needed to produce the streamlines

load("outputs/main_output.mat") % Load in the outputs and parameters from the main code

%% Remove the non-physical flows

% Get the coordinates of the appendage centers
pcenter1 = [PRAx  + (1+dsep/2)*cos(psi), PRAy];
pcenter2 = [PRAy  - (1+dsep/2)*cos(psi), PRAy];
pcenter3 = [-PRAx + (1+dsep/2)*cos(psi), PRAy];
pcenter4 = [-PRAx - (1+dsep/2)*cos(psi), PRAy];

% Preallocate space for a mask.
mask = ones(size(Uflowx)); % Note: this assumes length(x) = length(y).

% For each of the points in the solution space
for i = 1:length(x)
    for j = 1:length(y)

        % Get the current point
        p = [x(i),y(j)];

        %determine if the point is "in" an appendage
        a = sign(norm(p-pcenter1)-1);
        b = sign(norm(p-pcenter2)-1);
        c = sign(norm(p-pcenter3)-1);
        d = sign(norm(p-pcenter4)-1);

        % If the point is too close to an appendage remove it
        if a+b+c+d<4
            mask(i,j)=0;
        end
    end
end

% Remove all the non-physical points
Ux = Uflowx.*mask;
Uy = Uflowy.*mask;

%Umag = sqrt(Ux.^2 + Uy.^2); % Get flow field magnitude
%imagesc(x,y,Umag)

%% Simulate a particle which is subject to the remaining flows
% % Forward Euler
% 
% T = 50; % Total time.
% dt = 0.1; % Timestep.
% N = T/dt; % Number of timesteps.
% 
% scatter(stks(:,1),stks(:,2)); hold on
% 
% for x0 = linspace(5,6.5,100)
% 
% pos = zeros(N,2);
% pos(1,:) = [x0,12.5]; % particle initial position
% 
% for t = 2:N
% 
%     [~,indx] = min(abs(pos(t-1,1)-x));
%     [~,indy] = min(abs(pos(t-1,2)-y));
% 
%     pos(t,:) = pos(t-1,:) + [Uflowx(indx,indy),Uflowy(indx,indy)]*dt;
% 
% end
% 
% plot(pos(:,1),pos(:,2))
% hold on
% 
% end

%% Simulate a particle which is subject to the remaining flows
% ODE45 - Runge-Kutta 45

% Parameters for the solver
T = 100; % Total time.
dt = 0.1; % Timestep.
N = T/dt; % Number of timesteps.
times = linspace(0,T,N); % Time points at which to solve.
Npts = 100; % Number of particles to solve streamlines for.
rad = 1.1; % Radius of the particle for which the streamline is being solved
c = jet(Npts); % Colormap for each particle.
eps = 10; % Interaction strength scaling
imethod = 'spline'; % Set the interpolation calculation method to be used.

% Get the points in the pre-solved space on a grid
[gridx,gridy] = ndgrid(x,y);

% Get the interpolators for each component of the flow
UX = griddedInterpolant(gridx,gridy,flipud(imrotate(Uflowx,90)),imethod);
UY = griddedInterpolant(gridx,gridy,flipud(imrotate(Uflowy,90)),imethod);

% Define the equation of motion for the 
funct = @(t,x) [UX(x(2),x(1));UY(x(2),x(1))] + [DUX(x(1),x(2),rad,eps,psi,PRAx,PRAy,dsep);DUY(x(1),x(2),rad,eps,psi,PRAx,PRAy,dsep)]; % Steric interactions
%funct = @(t,x) [UX(x(2),x(1));UY(x(2),x(1))]; % No steric interactions

x0 = linspace(-8,8,Npts); % Vary the initial condition

% Preallocate space for the trajectory of each of the particles
trajx = zeros(length(times),Npts);
trajy = zeros(length(times),Npts);

% Main loop
for i = 1:Npts % Loop through the particles

    ic = [x0(i);10]; % Get the initial conditions.
    [tout,out] = ode45(funct, times', ic);  % Solve for the trajectory.
    
    % Store the outcome of the streamine.
    trajx(:,i) = out(:,1);
    trajy(:,i) = out(:,2);

end

% axis equal

%% Plot the outputs of the trajectory calculations

close all % Clear the current figure space 

% Prepare the file space for storing the outputs
delete movies/*
rmdir movies
mkdir movies/

rad = 1.1;
theta = linspace(0,2*pi,100);

for n = 1:5:length(times)/2

    hold off;

    scatter(stks(:,1),stks(:,2),5); hold on

    for i = 1:length(trajx(n,:))
        plot(trajx(n,i) + rad*cos(theta),trajy(n,i)+ rad*sin(theta),'lineWidth',5,'color',c(i,:));
    end

    ylim([-12,15])
    axis equal

    saveas(gcf,['movies/' num2str(n) '.png'])

end

%%

% Specify the folder containing the PNG images
imageFolder = 'movies/';

% Create a VideoWriter object to define the output video
outputVideo = VideoWriter('output_movie.mp4', 'MPEG-4');

% Set the desired frame rate (adjust as needed)
outputVideo.FrameRate = 10;

% Open the video writer object
open(outputVideo);

% Loop through each image and write it to the video
for i = 1:5:length(times)/2
    % Read the current image
    thisImage = imread(['movies/' num2str(i) '.png']);
    
    % Write the image to the video
    writeVideo(outputVideo, thisImage);
end

% Close the video writer object to finish the file
close(outputVideo);

% Display a message
disp('Movie creation completed!');

%%

function [dux] = DUX(x,y,rad,eps,psi,PRAx,PRAy,dsep)

    psi = 0;
    PRAx = 5;
    PRAy = 5;
    dsep = 1;
    rad = 1.1;
    eps = 10;

    p = [x,y];
    pcenter1 = [PRAx  + (1+dsep/2)*cos(psi), PRAy];
    d1 = norm(p-pcenter1);
    pcenter2 = [PRAy  - (1+dsep/2)*cos(psi), PRAy];
    d2 = norm(p-pcenter2);
    pcenter3 = [-PRAx + (1+dsep/2)*cos(psi), PRAy];
    d3 = norm(p-pcenter3);
    pcenter4 = [-PRAx - (1+dsep/2)*cos(psi), PRAy];
    d4 = norm(p-pcenter4);

    rel = min([d1,d2,d3,d4]);

    if 1*(0.5+rad - rel) > 0

        switch rel

            case d1
                dux = eps*(0.5+rad - d1)*(x-pcenter1(1))/d1;
            case d2
                dux = eps*(0.5+rad - d2)*(x-pcenter2(1))/d2;
            case d3
                dux = eps*(0.5+rad - d3)*(x-pcenter3(1))/d3;
            case d4
                dux = eps*(0.5+rad - d4)*(x-pcenter4(1))/d4;
        end

    else
        dux = 0;
    end

end

function [duy] = DUY(x,y,rad,eps,psi,PRAx,PRAy,dsep)


    p = [x,y];
    pcenter1 = [PRAx  + (1+dsep/2)*cos(psi), PRAy];
    d1 = norm(p-pcenter1);
    pcenter2 = [PRAy  - (1+dsep/2)*cos(psi), PRAy];
    d2 = norm(p-pcenter2);
    pcenter3 = [-PRAx + (1+dsep/2)*cos(psi), PRAy];
    d3 = norm(p-pcenter3);
    pcenter4 = [-PRAx - (1+dsep/2)*cos(psi), PRAy];
    d4 = norm(p-pcenter4);

    rel = min([d1,d2,d3,d4]);

    if 1*(0.5+rad - rel) > 0

        switch rel
            case d1
                duy = eps*(0.5+rad - d1)*(y-pcenter1(2))/d1;
            case d2
                duy = eps*(0.5+rad - d2)*(y-pcenter2(2))/d2;
            case d3
                duy = eps*(0.5+rad - d3)*(y-pcenter3(2))/d3;
            case d4
                duy = eps*(0.5+rad - d4)*(y-pcenter4(2))/d4;
        end

    else
        duy = 0;
    end

end
