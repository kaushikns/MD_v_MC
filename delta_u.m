function [du] = delta_u(r,r0,ind,L,dim)
n = length(r);
u1 = 0;
u2 = 0;
for i=1:n
    if i~=ind
        rij = r(i,:) - r(ind,:);
        for k=1:dim
            if rij(k) > 0.5*L
                rij(k) = rij(k)-L;
            elseif rij(k) < -0.5*L
                rij(k) = rij(k)+L;
            end
        end
        mag_rij = sqrt(sum(rij.^2));
        u1 = u1 + 4*((1/mag_rij)^12 - (1/mag_rij)^6);
        
        rij = r(i,:) - r0(ind,:);
        for k=1:dim
            if rij(k) > 0.5*L
                rij(k) = rij(k)-L;
            elseif rij(k) < -0.5*L
                rij(k) = rij(k)+L;
            end
        end
        mag_rij = sqrt(sum(rij.^2));
        u2 = u2 + 4*((1/mag_rij)^12 - (1/mag_rij)^6);
    end
end
du = u1 - u2;
