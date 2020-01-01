function [g,r0] = gr(r,n,L,dim)
dens = n/L^dim;
dr = L/100;
r0 = (0:dr:(0.5*L))';
g = zeros(length(r0)-1,1);
m = n*(n-1);
rij = zeros(m,1);
m = 1;
for i=1:n
    for j=1:n
        if i~=j
        delta = r(i,:) - r(j,:);
         for k=1:dim
            if delta(k) > 0.5*L
                delta(k) = delta(k)-L;
            elseif delta(k) < -0.5*L
                delta(k) = delta(k)+L;
            end
         end
        %delta = mod(delta,L);
        rij(m) = sqrt(dot(delta,delta));
        m = m + 1;
        end
    end
end
m = m-1;
for i=2:length(r0)
    count = 0;
    for j=1:m 
        if rij(j)>=r0(i-1) && rij(j) < r0(i)
            count = count + 1;
        end
    end
    if dim == 1
        dens_act = count/(n*(r0(i)^dim - r0(i-1)^dim));
    elseif dim == 2
        dens_act = count/(n*pi*((r0(i))^dim - (r0(i-1))^dim));
    elseif dim == 3
        dens_act = 3*count/(4*n*pi*((r0(i))^dim - (r0(i-1))^dim));
    end
    g(i) = dens_act/dens;
end