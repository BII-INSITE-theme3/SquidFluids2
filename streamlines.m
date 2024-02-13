% Title: Streamline calculations.
% Author: Stephen Williams.
% Notes: 1. Code is non-dimensionalised. 
%--------------------------------------------%

clear all
close all

%% 

save("U_streams","Uflowx","Uflowy","x","y","stks")

%% 

load("U_streams.mat")
parameters

%%

pcenter1 = [PRAx  + (1+dsep/2)*cos(psi), PRAy];
pcenter2 = [PRAy  - (1+dsep/2)*cos(psi), PRAy];
pcenter3 = [-PRAx + (1+dsep/2)*cos(psi), PRAy];
pcenter4 = [-PRAx - (1+dsep/2)*cos(psi), PRAy];

mask = ones(size(Uflowx));

for i = 1:length(x)
    for j = 1:length(y)

        p = [x(i),y(j)];

        a = sign(norm(p-pcenter1)-1);
        b = sign(norm(p-pcenter2)-1);
        c = sign(norm(p-pcenter3)-1);
        d = sign(norm(p-pcenter4)-1);

        if a+b+c+d<4
            mask(i,j)=0;
        end
    end
end

Ux = Uflowx.*mask;
Uy = Uflowy.*mask;

Umag = sqrt(Ux.^2 + Uy.^2); % Get flow field magnitude
imagesc(x,y,Umag)

%%

T = 50; % Total time.
dt = 0.1; % Timestep.
N = T/dt; % Number of timesteps.

scatter(stks(:,1),stks(:,2)); hold on

for x0 = linspace(5,6.5,100)

pos = zeros(N,2);
pos(1,:) = [x0,12.5]; % particle initial position

for t = 2:N

    [~,indx] = min(abs(pos(t-1,1)-x));
    [~,indy] = min(abs(pos(t-1,2)-y));

    pos(t,:) = pos(t-1,:) + [Uflowx(indx,indy),Uflowy(indx,indy)]*dt;

end

plot(pos(:,1),pos(:,2))
hold on

end
