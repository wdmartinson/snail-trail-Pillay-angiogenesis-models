function [BCModel_StalkCellDensity, BCModel_TipCellDensity, Omega, TimeMesh] = SnailTrail_1D_PDE(TC_InitCond, EC_InitCond)
% Solves the 1D snail-trail on the [0, 1] domain using the Method
% of Lines:
% n_t = D*n_xx-chi*(n*c_x)_x + lambda*n*c- beta_e*n*e-beta_n*n.^2
% e_t = kappa/h*abs(D*n_x - chi*c_x*n)
% 
% Parameter Values in Paper:
% chi = 0.4; beta_n = 160; D = 1e-3; lambda = 0.16; c_x = 1; beta_e =
% 5.0648.
% Initial conditions are ABM distributions from CA model, averaged in the
% y-direction, at t = 0.2 (and are linearly interpolated onto the PDE mesh
% Omega).
%--------------------------------------------------------------------------
%% Set up Domain, Initialize Parameters, Set Initial Condition
% Space-Time Domain:
Omega = linspace(0, 1, 201)'; % X-Coordinates
TimeMesh = 0:1/160:2;
N = length(Omega);
Dx = Omega(2);

% Parameters:
k = 100;
h = 1/200;
beta_n = 160; chi = beta2*k*h^2; D = beta2*h^2/4; beta_e = 1e-12*beta2; lambda = 0.16; kappa = 2;

% Initial Condition, specified as a vector
z0 = zeros(2*N, 1); % [p; n; new_p determined by mu*N approximation; new_n determined by mu*N approximation for p]

% CA Model IC at t = eps, where eps > 0
z0(1:N) = EC_InitCond;
z0(N+1:2*N) = TC_InitCond;

% TAF Field
c = linspace(0,max(Omega),N)'; % c(x) = x

e = ones(N,1);

% Use Block Sparse Matrices to define matrix system for MOL ode below:
N_N = spdiags([(D/Dx/Dx+chi/2/Dx)*e, (-2*D/Dx/Dx)*e, (D/Dx/Dx-chi/2/Dx)*e], -1:1, N, N);

% Establish BCs
% Neumann BC @ x = 0
N_N(1,1) = -2*D/Dx/Dx-2*chi/Dx-chi^2/D;
N_N(1,2) = 2*D/Dx/Dx;

% Neumann BC @ x = 1
N_N(end,end-1) = 2*D/Dx/Dx;
N_N(end,end) = -2*D/Dx/Dx + 2*chi/Dx-chi^2/D;

E_N = spdiags([(-D/Dx-chi)*e, (D/Dx)*e], 0:1, N, N);
E_N(1,1) = 0; % No flux at boundary
E_N(end,end) = 0; % No flux at boundary

A = [sparse(N,N), E_N;...
    sparse(N,N), N_N];
%% Solve the PDE using the Method of Lines
opts = odeset('MaxStep', min([Dx/chi,Dx^2/2/D]));
[~, Sols] = ode15s(@MOL_ODE_Byrne_Chaplain_Model, TimeMesh, z0, opts);
BCModel_StalkCellDensity = Sols(:, 1:N)';
BCModel_TipCellDensity = Sols(:, N+1:2*N)';
%--------------------------------------------------------------------------
% Subfunctions
    function dzdt = MOL_ODE_Byrne_Chaplain_Model(~, z)
        % Uses central differencing for n_x terms, central
        % differencing for 2nd order derivatives
        p = z(1:N); n = z(N+1:2*N); 
        b = zeros(2*N,1);
        b(N+1:2*N) = lambda*n.*c - beta_e*n.*p - beta_n*n.^2;
        
        dzdt = A*z + b;
        dzdt(1:N) = kappa/h*abs(dzdt(1:N));
    end % function MOL_ODE_Byrne_Chaplain_Model
end % function Byrne_Chaplain_Model