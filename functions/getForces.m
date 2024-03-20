% Function to solve the inverse problem to get the stokeslet forces.

function [F] = getForces(stks1,eps1,U01)

    % Create the array of boundary velocities (from BVP)
    [nStok,~] = size(stks1); % Get number of stokeslets required at given linear density.
    BdryVelo = zeros(nStok,2); % Preallocate an array to store the boundary velocities (for BVP).

    %% No-slip boundaries -- stks(:,3) == 1, this code currently does nothing, so is commented out.
    %ind = find(stks(:,3)==1);
    %BdryVelo(ind,:) = 0; % Set zero-flow on the channel boundaries

    %% Poisuelle boundary flow -- stks(:,3) == 2
    % Flow going as 1-(r/a)^2, r = distance from channel center, a = half channel-width.
    ind = find(stks1(:,3)==2); % Get relevant Stokeslets for which this boundary condition holds
    nTemp = length(ind); % Get the number of Stokeslets contained in this set.
    BdryVelo(ind,:) = poisuelleFlow(nTemp,U01); % Prescibe the boundary velocity at these points.

    %% No-slip boundaries -- stks(:,3) == 3
    % Can add other boundary conditions for the bottom of a channel here.

    %% No-slip boundaries -- stks(:,3) == 4,5,6,7
    % Boundary conditions following those found in Nawroth 2017, see SM.

    ind = find(stks1(:,3)==4); nTemp = length(ind);
    BdryVelo(ind,:) = surfaceFlow(nTemp,floor(0.75*nTemp));

    ind = find(stks1(:,3)==5); nTemp = length(ind);
    BdryVelo(ind,:) = surfaceFlow(nTemp,floor(0.5*nTemp));

    ind = find(stks1(:,3)==6); nTemp = length(ind);
    BdryVelo(ind,:) = surfaceFlow(nTemp,floor(0.5*nTemp));

    ind = find(stks1(:,3)==7); nTemp = length(ind);
    BdryVelo(ind,:) = surfaceFlow(nTemp,floor(0.75*nTemp));

    %% Create the linear system to solve, part 1.
    % Put the velocity boundary conditions in vertical form.

    BdryVertical = zeros(nStok*2,1);
    BdryVertical(1:2:end) = BdryVelo(:,1); BdryVertical(2:2:end) = BdryVelo(:,2);
    % quiver(stks1(:,1),stks1(:,2),BdryVertical(1:2:end),BdryVertical(2:2:end)); % Optional: View the constructed flow profile

    %% Create the linear system to solve, part 2.
    % Create a matrix corresponding to the stokeslet for the full system.

    S = zeros(2*nStok); % Preallocate the full linear stokeslet.

    % Find the forces.
    for l = 1:2*nStok % Loop through the stokeslet-components
        for k = 1:2*nStok % Loop through the stokeslet-components
            %
            n1 = ceil(k/2); % Get the stokelets number of the influence stokeslet.
            n2 = ceil(l/2); % Get the stokelets number of the influenced stokeslet.
            %
            tempStks1 = stks1(n1,1:2); % Get the position of the influence stokeslet.
            tempStks2 = stks1(n2,1:2); % Get the position of the influenced stokeslet.
            %
            r = sqrt(norm(tempStks1-tempStks2)^2 + eps1^2) + eps1; % Get the reg `distance' between them.
            rho = (r+eps1)/(r*(r-eps1)); % Get the "rho", for convenience.
            %
            S(l,k) = -(log(r)-eps1*rho)*(mod(k,2)==mod(l,2)) ... % Get the log term contribution.
                + (tempStks1(2-mod(k,2))-tempStks2(2-mod(k,2)))*(tempStks1(2-mod(l,2))-tempStks2(2-mod(l,2)))*rho/r; % Get the <....> contribution.
        end
    end

    %% Create the linear system to solve, part 2.
    % Solve for the forces required to satisfy the system.

    % Calculate the forces and then put them into vector format for output.
    % ForceVertical = zeros(2*nStok,1);
    %tol = 1e-4;
    %maxit = 100;
    %[ForceVertical,~] = gmres(S,BdryVertical,[],tol,maxit); % Inv here is inefficient, but seems to work well. It is likely poorly conditioned for very large numbers of Stokelets.
    ForceVertical = S\BdryVertical; % Backslash method, potential source of error.
    F = zeros(nStok,2);
    F(:,1) = ForceVertical(1:2:end);
    F(:,2) = ForceVertical(2:2:end);
    
    % Set the net force to zero
    F(:,1) = F(:,1) - mean(F(:,1));
    F(:,2) = F(:,2) - mean(F(:,2));

    % Visualise the forces calculated. May need to remove large values to aid the viewing.
    % quiver(stks1(:,1),stks1(:,2),F(:,1),F(:,2))

end

