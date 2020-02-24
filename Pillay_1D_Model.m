function [PillayModel_StalkCellDensity, PillayModel_TipCellDensity, Omega, TimeMesh] = Pillay_1D_Model(TC_InitCond, EC_InitCond)
% Solves the 1D Pillay PDE using the Method of Lines
% 
% n_t = (1-a_n*n-a_e*e)*D*n_xx-chi*(n*c_x)_x + lambda*n*c- a_e*mu*n*e-a_n*mu*n.^2 
% e_t = mu*n + mu*n^2 + a_n*n*(D*n_xx-chi*(n*c_x)_x)
%
% (Note: Equation for c not used because c(x) is assumed linear)
% 
% Parameter Values in Original Paper:
% chi = 0.4; mu = 160 ; D = 1e-3; lambda = 0.16; a_n = 1; a_e = 0.0343;
% Initial conditions are ABM distributions from CA model, averaged in the
% y-direction, at t = 0.2 (and are linearly interpolated onto the PDE mesh
% Omega).
%--------------------------------------------------------------------------
%% Set up Domain, Initialize Parameters, Set Initial Condition
% Space-Time Domain:
Omega = linspace(0, 1, 201)'; % X-Coordinates
TimeMesh = 0.2:1/160:2;
N = length(Omega);
Dx = Omega(2);
% Parameters:
k = 100;
h = 1/200;
mu = 160; chi = mu*k*h^2; D = mu*h^2/4; lambda = 0.16; a_n = 1; a_e = 0.0343;
% Initial Condition, specified as a large vector
z0 = zeros(2*N,1); %[p; n; c] = z0

% CA Model IC from t = eps, where eps > 0
z0(1:N) = EC_InitCond;
z0(N+1:end) = TC_InitCond;

c = linspace(0,max(Omega),N)'; % c(x) = x

e = ones(N,1);
% Use Block Sparse Matrices to define matrix system for MOL ode below:
A = spdiags([(D/Dx/Dx+chi/2/Dx)*e, (-2*D/Dx/Dx)*e, (D/Dx/Dx-chi/2/Dx)*e], -1:1, N, N);
A(1,1) = -2*D/Dx/Dx-2*chi/Dx-chi^2/D;
A(1,2) = 2*D/Dx/Dx;
A(end,end-1) = 2*D/Dx/Dx;
A(end,end) = -2*D/Dx/Dx + 2*chi/Dx-chi^2/D;

%% Solve the PDE using the Method of Lines
opts = odeset('MaxStep', min([Dx/chi,Dx^2/2/D]));
[~, Sols] = ode15s(@MOL_ODE_Pillay_Model, TimeMesh, z0, opts);
PillayModel_StalkCellDensity = Sols(:, 1:N)';
PillayModel_TipCellDensity = Sols(:, N+1:2*N)';
%--------------------------------------------------------------------------
%% Subfunctions
    function dzdt = MOL_ODE_Pillay_Model(~, z)
        % Uses upwinding in x-direction for c_x, n_x terms, central
        % differencing for 2nd order derivatives
        p = z(1:N); n = z(N+1:2*N);
        dndt = (1-a_n*n-a_e*p).*(A*n) + lambda*n.*c - mu*(a_n*n.^2+a_e*n.*p);
        dpdt = mu*n+a_n*n.*(mu*n+A*n);
        
        dzdt = [dpdt; dndt];
    end % function MOL_ODE_Byrne_Chaplain_Model
end % function Byrne_Chaplain_Model