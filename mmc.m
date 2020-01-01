function [r,g,rg,u_out,P,t_mc] = mmc(n,L,T,dim,nStep)
u_out = zeros(nStep,1);
P = zeros(nStep,1); % Used interchangeably for pressure (2-3d) and MSD (1d)
beta = 1/T;
tol = 100;

% LJ simulation always done in 3d for this project
% 2d or 3d sim. of LJ fluid
if dim > 1
    [r0,~] = init3d(n,L,T,tol);
    r = r0;
    delta = 0.1*L;
    u = lj(r,L,dim);
end

% MC as a generator of overdamped Lang. dynamics (always 1d)
% 1d sim, free diffusion, BHO
if dim == 1
    r0 = zeros(n,1);
%     r0 = 5*ones(n,1);
    r = r0;
    r_ini = r0;
    % Free diffusion
    % u = 0;
    % BHO
    u = spring(n,r);
    delta = 0.01;
end

count = 0;
tic;
for i=1:nStep
    %"Sweep"
    p_acc = 0;
    for j=1:n
        q = rand;
        ind = ceil(q*n);
        dr = 2*delta*(rand(1,dim)-0.5);
        r(ind,:) = r0(ind,:) + dr;
        
        if dim == 1
            % Free diffusion
            % du = 0;
            
            % Overdamped BHO
            du = 0.5*(r(ind)*r(ind) - r0(ind)*r0(ind)); %Spring k = 1
        end
        % Periodic BC
        if dim > 1
            r(ind,:) = mod(r(ind,:),L);
            du = delta_u(r,r0,ind,L,dim);
        end
        
        u = u + du;
        accept = min(1,exp(-beta*du));
        s = rand;
        if s > accept
            r(ind,:) = r0(ind,:);
            p_acc = p_acc-1;
            u = u - du;
        end
        
        % MSD computation for 1d case
        if dim == 1
            dx = r(ind,1) - r_ini(ind,1);
            dx2 = dx*dx;
            P(i) = P(i) + dx2;
        end
        r0 = r;
        p_acc = p_acc + 1;
    end
    
    g = 0;
    rg = g;
    
    if dim > 1
        %%% Modify move size if acceptance neq 0.5
        p_acc = p_acc/n;
        % disp(p_acc);
        if p_acc > 0.5
            delta = delta*1.1;
            if delta > L
                delta = L/1.1;
            end
        elseif p_acc < 0.5
            delta = delta*0.9;
        end
        %%%
        
        %%% Radial dist function
        if mod(i,10)==0 && i > 0.5*nStep
            count = count+1;
            [g,rg] = gr(r,n,L,dim);
            if count==1
                g_old = g;
                continue;
            end
            g = ((count-1)*g_old + g)/count;
            g_old = g;
            div = 3;
            dens_hist(r,n,L,div);
        end
        %%%
    end
    if dim == 1 && i == nStep/20
        r_mid = r;
    end
    u_out(i) = u;
    %     P(i) = pres(n,L,T,r);
end
if dim == 1
    P = P/n;
    r = [r_mid, r];
end
t_mc = toc;