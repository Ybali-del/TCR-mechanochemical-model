% -------------------------------------------------------------------------
% Script: saddle_point_analysis_energy_landscape.m
% Purpose: Identify local minima and saddle points in the force-loaded 
%          mechanochemical energy landscape V_f(L, d).
% Paper:   RSOS submission: "Mechanochemical Energy Landscapes Under Force 
%          Explain Catch–Slip Bonds in T-Cell Activation"
% -------------------------------------------------------------------------

%% ------------------------ CONSTANTS ---------------------------------
theta_0 = 0.60 * pi;      % radians
theta_1 = 0.80 * pi;      % radians
c       = 0;              % pN*nm
k_theta = 80.9;         % pN*nm
D_0     = 40.30;          % pN*nm
d_0     = 0.151;          % nm
sigma   = 0.8095* pi;    % radians

theta = @(L) 2 * asin(L / 22); 
V_theta = @(L) 0.5 * k_theta * (theta(L) - theta_0).^2; 
D_theta = @(L) D_0 * exp(-((theta(L) - theta_1).^2)/(2*sigma^2)) + c; 
B = @(L,d) D_theta(L).* ((1 - exp(-d/d_0)).^2 - 1);

for F = 1:5
    fprintf('\n================ FORCE F = %d =================\n', F);

    V_f = @(L,d) V_theta(L) + B(L,d) - F*(L + d);

    grad_Vf = @(x)[
        (V_f(x(1)+1e-6,x(2)) - V_f(x(1)-1e-6,x(2))) / (2e-6);
        (V_f(x(1),x(2)+1e-6) - V_f(x(1),x(2)-1e-6)) / (2e-6)
    ];

    critical_points = [];

    for i = 1:10
        L_init = 22 + (22 - 5)*rand();
        d_init = -0.4 + (1 + 0.3)*rand();

        [sol,~,flag] = fsolve(grad_Vf, [L_init d_init],optimoptions('fsolve','Display','none'));

        if flag > 0
            if isempty(critical_points) || all(vecnorm(critical_points - sol,2,2) > 1e-3)
                critical_points = [critical_points; sol];
            end
        end
    end

    fprintf('Critical points and Hessian classification:\n');

    for i = 1:size(critical_points,1)
        L = critical_points(i,1);
        d = critical_points(i,2);

        h=1e-6;
        d2L = (V_f(L+h,d)-2*V_f(L,d)+V_f(L-h,d))/h^2;
        d2d = (V_f(L,d+h)-2*V_f(L,d)+V_f(L,d-h))/h^2;
        dLd = (V_f(L+h,d+h)-V_f(L+h,d-h)-V_f(L-h,d+h)+V_f(L-h,d-h))/(4*h^2);

        H=[d2L dLd; dLd d2d];
        eigs_H=eig(H);

        if all(eigs_H>0)
            nature='Local Minimum';
        elseif all(eigs_H<0)
            nature='Local Maximum';
        else
            nature='Saddle Point';
        end

        fprintf('L=%.4f, d=%.4f  | %s\n', L, d, nature);
        fprintf('Eigenvalues: %.4f, %.4f\n', eigs_H(1), eigs_H(2));
    end
end
