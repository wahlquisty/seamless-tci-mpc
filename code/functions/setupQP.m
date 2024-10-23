function [H, fT, A, b] = setupQP_u(Phi, Gamma, alpha_u, umax, r, N, x0)
%SETUPQP_U Set up the Quadratic Programming (QP) problem for solving control input.
%
%   [H, fT, A, b] = SETUPQP_U(Phi, Gamma, alpha_u, umax, r, N, x0) 
%   sets up the matrices required for solving the QP problem . 
%
%   Inputs:
%   - Phi: Discrete state matrix.
%   - Gamma: Control input matrix.
%   - alpha_u: Weighting factor for the control input rate penalty.
%   - umax: Maximum control input.
%   - r: Reference trajectory.
%   - N: Prediction horizon (number of future steps).
%   - x0: Initial state vector.
%
%   Outputs:
%   - H: Hessian matrix for the QP.
%   - fT: Linear term in the QP objective function.
%   - A: Matrix for the inequality constraints (control input limits).
%   - b: Vector for the inequality constraints (upper/lower control limits).

    n = length(x0);
    Phij = eye(n); % Phi^j
    Phij1Gamma = zeros(N,1); % Row j will hold Phi^j_1*Gamma
    Phij4Gamma = zeros(N,1); % Row j will hold Phi^j_4*Gamma
    E1 = zeros(N, n);
    E4 = zeros(N, n);
    F1 = zeros(N, N);
    F4 = zeros(N, N);

    for j = 1:N
        Phij1Gamma(j) = Phij(1,:) * Gamma;
        Phij4Gamma(j) = Phij(4,:) * Gamma;
        for i = 1:j
            F1(j,i) = Phij1Gamma(j-i+1);
            F4(j,i) = Phij4Gamma(j-i+1);
        end
        Phij = Phij * Phi;
        E1(j, :) = Phij(1,:);
        E4(j, :) = Phij(4,:);
    end

    % For constraint on u(k)-u(k-1)
    V = zeros(N,N);
    for j = 1:N-1
        V(j,j) = -1;
        V(j,j+1) = 1;
    end
    V(j,j) = -1;

    D = eye(N);

    % D(end, end) = alpha
    % Speed-up is possible here:
    % 1. H is symmetric, so only need to compute half
    % 2. D and F1 are both sparse
    DF4 = D * F4;                    % Compute only once
    aVtV = alpha_u*(V')*V;

    % H=F1'*DF1                  % Compute only once
    H = sparse(0.5 * (F4' * DF4 + DF4' * F4 + aVtV));     % Same but numerically more robust
    % H = (H + H') / 2; % Ensure H is symmetric % ADDED

    f0T = E4' * DF4;                 % Compute only once
    f1T = r' * DF4;                  % Recompute in every iteration
    fT = x0' * f0T - f1T;              % Recompute in every iteration

    A = sparse([-eye(N); eye(N)]);            % Compute only once
    b = [zeros(N, 1); umax*ones(N, 1)]; % Recompute lower part in each iteration
end