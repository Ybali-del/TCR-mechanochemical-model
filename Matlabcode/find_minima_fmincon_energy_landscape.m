% -------------------------------------------------------------------------
% Script: find_minima_fmincon_energy_landscape.m
% Purpose: Locate local minima of the mechanochemical energy landscape 
%          V_f(L, d) using fmincon for forces F = 1:5.
% -------------------------------------------------------------------------

for F = 1:5
    fprintf('\n================ FORCE F = %d =================\n', F);
 theta_0 = 0.60 * pi;      % radians
theta_1 = 0.80 * pi;      % radians
c       = 0;              % pN*nm
k_theta = 80.9;         % pN*nm
D_0     = 40.30;          % pN*nm
d_0     = 0.151;          % nm
sigma   = 0.8095* pi;    % radians

    theta = @(L) 2*asin(L/22);
    V_theta=@(L) 0.5*k_theta*(theta(L)-theta_0).^2;
    D_theta=@(L) D_0*exp(-((theta(L)-theta_1).^2)/(2*sigma^2)) + c;
    B=@(L,d) D_theta(L).* ((1 - exp(-d/d_0)).^2 - 1);
    V_f=@(x) V_theta(x(1))+B(x(1),x(2)) - F*(x(1)+x(2));

    lb=[0 0]; ub=[22 2];

    seeds = [0.1 0.1; 8 0.3; 10 0.1; 22 2];
    results=[];

    for k=1:size(seeds,1)
        [x_opt, fval] = fmincon(V_f,seeds(k,:),[],[],[],[],lb,ub);
        results=[results; x_opt fval];
    end

    results=unique(round(results,6),'rows');
    disp(results);

    for i=1:size(results,1)
        H = compute_hessian(V_f, results(i,1:2));
        eigs_H = eig(H);

        if all(eigs_H>0)
            disp('Local Minimum Found:')
        end
        disp(results(i,:))
    end
end

function H = compute_hessian(f,x)
    eps=1e-5;
    n=length(x);
    H=zeros(n);

    for i=1:n
        for j=1:n
            x_ij=x; x_ij(i)=x_ij(i)+eps; x_ij(j)=x_ij(j)+eps;
            x_i=x; x_i(i)=x_i(i)+eps;
            x_j=x; x_j(j)=x_j(j)+eps;

            H(i,j) = (f(x_ij)-f(x_i)-f(x_j)+f(x))/eps^2;
        end
    end
end
