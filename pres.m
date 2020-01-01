function [P] = pres(n,L,T,r)
V = L^3;
P = n*T;
for i=1:n
    for j=i+1:n
        rij = r(i,:) - r(j,:);
            for k=1:3
                if rij(k) > 0.5*L
                    rij(k) = rij(k)-L;
                elseif rij(k) < -0.5*L
                    rij(k) = rij(k)+L;
                end
            end
            mag_rij = sqrt(sum(rij.^2));
            fij = 24*(2/mag_rij^14 - 1/mag_rij^8)*rij;
            P = P + dot(rij,fij)/3;
    end
end
P = P/V;