function [f] = force(r,L,dim)
n = length(r);
f = zeros(n,dim);
% 1d sim, BHO
if dim == 1
    f = -r; % Spring k = 1
    return;
end
% Continue with a 2d or 3d LJ simulation if dim neq 1
for i=1:n
    for j=1:n
        if i~=j
            rij = r(i,:) - r(j,:);
            for k=1:dim
                if rij(k) > 0.5*L
                    rij(k) = rij(k)-L;
                elseif rij(k) < -0.5*L
                    rij(k) = rij(k)+L;
                end
            end
            mag_rij = sqrt(sum(rij.^2));
            fij = 24*(2*mag_rij^(-14) - mag_rij^(-8))*rij;
            f(i,:) = f(i,:) + fij;
        end
    end
end
end

