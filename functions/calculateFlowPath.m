function [Uflowx,Uflowy] = calculateFlowPath(stks1,F1,x1,y1,eps1)

    Uflowx = zeros(length(x1),1);
    Uflowy = zeros(length(y1),1);
    [nStok,~] = size(stks1);
    
    for position = 1:length(x1)
    
        Stemp = zeros(2,2);
        tempStks = stks1(:,1:2);
        tempF = F1;
        tempX = x1;
        tempY = y1;
    
        p = [tempX(position),tempY(position)]'; % Get the position of consideration, i and j done make sense to me here??
    
        for n = 1:nStok
    
            pN = tempStks(n,:)'; % Get the position of stokeslet N.
            Ftemp = tempF(n,:)'; % Get the forces of stokeslet N.
            r = sqrt(norm(p - pN).^2 + eps1^2) + eps1; % Distance, considered to stokeslet N.
            rho = (r+eps1)/(r*(r-eps1)); % Rho, considered to stokeslet N.
    
            for k = 1:2
                for l = 1:2
                    Stemp(k,l) = -(log(r)-eps1*rho)*(k==l) + (p(k)-pN(k))*(p(l)-pN(l))*rho/r;
                end
            end
    
            U = Stemp*Ftemp;
            Uflowx(position) = Uflowx(position) + U(1);
            Uflowy(position) = Uflowy(position) + U(2);
    
        end
    
    end

end

