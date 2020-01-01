function [dx,r,v,Tk,g,rg,u_out,P,t_md] = md(n,L,dt,T,dim,nStep)
tol = 10^-6;
% LJ simulation always done in 3d for this project
if dim > 1
    [r,v] = init3d(n,L,T,tol);
    dx = zeros(nStep,n);
    
    % 1d sim, BHO
elseif dim == 1
    r = zeros(n,1);
    v = rand(n,1);
    v = v - mean(v);
end
r0 = r;
r_ini = r;
f = force(r,L,dim);

% Nose Hoover thermostat variables
p = 5;
Q = zeros(p,1);
tau = 25*dt;
Q(1) = n*T*tau^2;
for j = 2:p
    Q(j) = 3*T;
end

xi = zeros(p,1);
%Data analysis part
t = 0;
count = 0;
g = 0;
rg = 0;
iter = 0;
u_out = zeros(nStep,1);
Tk = sum(dot(v,v))/(dim*n);

% Used interchangeably for pressure (2-3d) and MSD (1d)
P = zeros(nStep,1);

tic;

% MSD computation for 1d case
if dim == 1
    for i=1:nStep
        dx = r - r_ini;
        P(i) = dot(dx,dx);
        f = force(r,L,dim);
        
        % MD
        r = r + dt*v;
        v = v + dt*f - dt*v + sqrt(2*dt)*normrnd(0,1,[n 1]);
        
        % BD
%          r = r + dt*f + sqrt(2*dt)*normrnd(0,1,[n 1]);
    end
    P = P/n;
    t_md = toc;
    return;
end

% Use only if you want to write config to file
% fileid = fopen('md.xyz','a+');

% Use only if you want to write U to file
% fileid = fopen('potential.out','a+');

for i=1:nStep
    r = r + dt*v + 0.5*dt*dt*f;
    
    % Diffusivity calc for PBC in 2d-3d case
    r = mod(r,L);
    
    %%% Data generation for diff. calc.
    delta = r(:,1) - r0(:,1);
    for j=1:n
        if delta(j) > 0.5*L
            delta(j) = delta(j) - L;
        elseif delta(j) < -0.5*L
            delta(j) = delta(j) + L;
        end
        if i==1
            dx(i,j) = delta(j);
        else
            dx(i,j) = dx(i-1,j) + delta(j);
        end
    end
    %%%
    
    
    %%% Radial dist function
    if mod(iter,10)==0 && iter > 0.5*nStep
        count = count+1;
        [g,rg] = gr(r,n,L,dim);
        if count==1
            g_old = g;
            continue;
        end
        g = ((count-1)*g_old + g)/count;
        g_old = g;
    end
    %%%
    
    r0 = r;
    
    v = v + 0.5*dt*f - 0.5*xi(1)*dt*v;
    f = force(r,L,dim);
    v = v + 0.5*dt*f - 0.5*xi(1)*dt*v;
    
    Tk = sum(dot(v,v))/(dim*n);
    xi_dot = nh(T,Q,xi,p,Tk,n,dim);
    xi = xi + dt*xi_dot;
    
    % disp(Tk);
    
    t = t + dt;
    iter = iter + 1;
%     P(iter) = pres(n,L,Tk,r);
    u = lj(r,L,dim);
    u_out(iter) = u;
%     fprintf(fileid,'%f \n',u);
    
    % Use only if you want to write config to file
%     sig = 3.4; % In Angstrom
%     fprintf(fileid, '%d \n \n',n);
%     for j=1:n
%         fprintf(fileid,'Ar %f %f %f \n',r(j,1)*sig,r(j,2)*sig,r(j,3)*sig);
%     end
end
% fclose(fileid);
t_md = toc;