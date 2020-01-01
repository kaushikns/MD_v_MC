function [xi,xi_dot] = nh(T,Q,xi,p,Tk,n,dim)
xi_dot = zeros(p,1);
xi_dot(1) = -xi(1)*xi(2);
xi_dot(1) = xi_dot(1) + (dim*n*Tk - dim*n*T)/Q(1);
for j=2:p-1
    xi_dot(j) = (Q(j-1)*xi(j-1)*xi(j-1) - T)/Q(j) - xi(j)*xi(j+1);
end
xi_dot(p) = (Q(p-1)*xi(p-1)*xi(p-1) - T)/Q(p);
end