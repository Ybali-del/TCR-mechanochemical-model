function out = veff_min_and_saddles(F)
% -------------------------------------------------------------------------
% Computes the minimum and saddle point of V_f(L, d) for force F.
% -------------------------------------------------------------------------

if nargin < 1
    F = 100;
end

 theta_0 = 0.60 * pi;      % radians
theta_1 = 0.80 * pi;      % radians
c       = 0;              % pN*nm
k_theta = 80.9;         % pN*nm
D_0     = 40.30;          % pN*nm
d_0     = 0.151;          % nm
sigma   = 0.8095* pi;    % radians


Lmin=0; Lmax=14;
dmin=0; dmax=0.3;
h=1e-6;

theta=@(L) 2*asin(max(min(L,14),0)/14);
V_theta=@(L) 0.5*k_theta*(theta(L)-theta_0).^2;
D_theta=@(L) D_0*exp(-((theta(L)-theta_1).^2)/(2*sigma^2)) + c;
B=@(L,d) D_theta(L).*((1 - exp(-d/d_0)).^2 - 1);
Vf=@(L,d) V_theta(L)+B(L,d)-F*(L+d);

gradV=@(x)[
    (Vf(x(1)+h,x(2)) - Vf(x(1)-h,x(2))) / (2*h);
    (Vf(x(1),x(2)+h) - Vf(x(1),x(2)-h)) / (2*h)
];

hessV=@(x)[
    (Vf(x(1)+h,x(2)) - 2*Vf(x(1),x(2)) + Vf(x(1)-h,x(2))) / h^2, ...
    (Vf(x(1)+h,x(2)+h) - Vf(x(1)+h,x(2)-h) - Vf(x(1)-h,x(2)+h) + Vf(x(1)-h,x(2)-h)) / (4*h^2); ...
    (Vf(x(1)+h,x(2)+h) - Vf(x(1)+h,x(2)-h) - Vf(x(1)-h,x(2)+h) + Vf(x(1)-h,x(2)-h)) / (4*h^2), ...
    (Vf(x(1),x(2)+h) - 2*Vf(x(1),x(2)) + Vf(x(1),x(2)-h)) / h^2 ...
];

%% ---------------- LOCAL MINIMA (fmincon) ----------------
Lgrid=linspace(Lmin,Lmax,5);
dgrid=linspace(dmin,dmax,5);
seeds=[];

for i=1:5
    for j=1:5
        seeds(end+1,:)=[Lgrid(i), dgrid(j)];
    end
end

opts=optimoptions('fmincon','Display','none');

mins=[];
for k=1:size(seeds,1)
    [x, fval, flag]=fmincon(@(z)Vf(z(1),z(2)), seeds(k,:), [],[],[],[], ...
        [Lmin dmin],[Lmax dmax],[],opts);
    if flag>0
        mins=[mins; x fval];
    end
end

mins=unique(round(mins,6),'rows');

minima=[];
for i=1:size(mins,1)
    H=hessV(mins(i,1:2));
    ev=eig(H);
    if all(ev>0)
        minima=[minima; mins(i,1:2), Vf(mins(i,1), mins(i,2))];
    end
end

%% ---------------- SADDLES (fsolve) ----------------
seedsB=seeds;
opts2=optimoptions('fsolve','Display','none');

crit=[];
for k=1:size(seedsB,1)
    [sol,~,flag]=fsolve(@(z)gradV(z),seedsB(k,:),opts2);
    if flag>0
        crit=[crit; sol(:)'];
    end
end

crit=unique(round(crit,6),'rows');

saddles=[];
for i=1:size(crit,1)
    H=hessV(crit(i,:));
    ev=eig(H);
    if any(ev>0)&&any(ev<0)
        saddles=[saddles; crit(i,:), Vf(crit(i,1), crit(i,2))];
    end
end

out.minima=minima;
out.saddles=saddles;
end
