function [r,v] = init3d(n,L,T,tol)
% fileid = fopen('md.xyz','a+');
r = zeros(n,3);
v = zeros(n,3);
n13 = n^(1/3);
count = 1;
for i=1:ceil(n13)
    for j=1:ceil(n13)
        for k=1:ceil(n13)
            r(count,1) = (i-0.5)*L/ceil(n13);
            r(count,2) = (j-0.5)*L/ceil(n13);
            r(count,3) = (k-0.5)*L/ceil(n13);
            v(count,1) = rand;
            v(count,2) = rand;
            v(count,3) = rand;
            count = count + 1;
        end
    end
end

r = r(1:n,:);
% r = rand(n,3)*L;
v = v(1:n,:);
v = v - mean(v);
Tk = sum(dot(v,v))/(3*n);
v = v*sqrt(T/Tk);
f = force(r,L,3);
dt = 0.0005;
tau = 2*dt;
i = 1;

% MMC does not require going close to set temp.
if tol >= 1
    return;
end

% Use only if you want to write config to file
% fprintf(fileid, '%d \n \n',n);
% sig = 3.4; %In Angstrom
% for j=1:n
%     fprintf(fileid,'Ar %f %f %f \n',r(j,1)*sig,r(j,2)*sig,r(j,3)*sig);
% end

% Use only if you want to write U to file
% fileid = fopen('potential.out','a+');

while i==1 || abs(Tk-T) > tol
    
%     u = lj(r,L,3);
    % Use only if you want to write U to file
%     fprintf(fileid,'%f \n',u);
    
    r = r + dt*v + 0.5*dt*dt*f;
    r = mod(r,L);
    v = v + 0.5*dt*f;
    f = force(r,L,3);
    v = v + 0.5*dt*f;
    Tk = sum(dot(v,v))/(3*n);
    lambda2 = 1 + dt*(T-Tk)/(Tk*tau);
    v = sqrt(lambda2)*v;
    i = i+1;
    % Use only if you want to write config to file
%     fprintf(fileid, '%d \n \n',n);
%     for j=1:n
%         fprintf(fileid,'Ar %f %f %f \n',r(j,1)*sig,r(j,2)*sig,r(j,3)*sig);
%     end
end
% fclose(fileid);