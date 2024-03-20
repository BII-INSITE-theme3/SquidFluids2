function [Uflow] = surfaceFlow(nTemp1,rot1)

    Uflow = zeros(nTemp1,2); % Preallocate the flow vector

    % Get the boundary flow magnitude
    magnitude = zeros(nTemp1,1); % Store for the magnitudes, which is a combined sinosidal and zero profile
    indRange = floor(0.75*nTemp1); % Get the number of points for which the sinosidal curve is present (3/4 of the surface)

    theta = linspace(0,2*pi,indRange+1); % Parameterise this surface to have a full sin curve
    theta = theta(1:end-1); % Remove the repeating point

    magnitude(1:indRange) = sin(theta); % Get the magnitude of the sin across this curve
    magnitude = rotateArray(magnitude,rot1); % "Rotate" this array to align with the desired location (see rotateArray.m for details).

    % Construct the tangent vector along the curve
    theta = linspace(0,2*pi,nTemp1+1); % Parameterise
    theta = theta(1:end-1); % Remove the repeat point

    Uflow = magnitude.*[-cos(theta'),sin(theta')]; % Get the resulting flow combining the tangent and the magnitude

    % Optional: Remove the average so that the net flow is zero.
    % Uflow(:,1) = Uflow(:,1) - mean(Uflow(:,1));
    % Uflow(:,2) = Uflow(:,2) - mean(Uflow(:,2));

end

