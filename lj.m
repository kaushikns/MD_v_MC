function [u] = lj(r,L,dim)
u = 0;
n = length(r);
for i=1:n
    for j=i+1:n
        rij = r(i,:) - r(j,:);
            for k=1:dim
                if rij(k) > 0.5*L
                    rij(k) = rij(k)-L;
                elseif rij(k) < -0.5*L
                    rij(k) = rij(k)+L;
                end
            end
        mag_rij = sqrt(sum(rij.^2));
        u = u + 4*((1/mag_rij)^12 - (1/mag_rij)^6);
    end
end