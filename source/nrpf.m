function [V, delta, Pslack, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
% Deniz Temurcu 261089503
% This function performs Full Newton–Raphson power flow (polar form)

% Our inputs:
% Y      is the bus admittance matrix (n x n, complex)
% is     is the slack bus index (scalar, 1-based)
% ipq    is the vector of PQ-bus indices (column)
% ipv    is the vector of PV-bus indices (column)
% Pg     is the vector of active power generation (MW, n x 1)
% Qg     is the vector of reactive power generation (Mvar, n x 1)
% Pd     is the vector of active power demand (MW, n x 1)
% Qd     is the vector of reactive power demand (Mvar, n x 1)
% V0     is the vector of initial/specified voltage magnitudes (p.u., n x 1)
% Sbase  is the MVA base (scalar)
% toler  is the convergence tolerance on mismatch (∞-norm, p.u.)
% maxiter is the maximum number of NR iterations (scalar)

% Our outputs:
% V      is the vector of bus voltage magnitudes (p.u., n x 1)
% delta  is the vector of bus voltage angles (rad, n x 1) with delta(is)=0
% Psl    is the slack-bus active generation (MW)
% Qgv    is the vector of PV-bus reactive generations (Mvar), same order as ipv
% N      is the number of iterations performed (scalar)
% time   is the CPU time (s) of the NR loop (scalar)
%

% make sure the sizes match up
n = size(Y,1); % get the number of buses
if size(Y,2) ~= n
    error('Y must be square.');
end
if any([numel(Pg),numel(Qg),numel(Pd),numel(Qd),numel(V0)] ~= n)
    error('Pg,Qg,Pd,Qd,V0 must all be length n.');
end
ipq = ipq(:); ipv = ipv(:); % ensure column vectors
if any(ipq == is) || any(ipv == is)
    error('ipq and ipv must not include the slack bus.');
end
if ~isempty(intersect(ipq, ipv))
    error('ipq and ipv must be disjoint.');
end

% initialize
Pinjection = (Pg(:) - Pd(:)) / Sbase; % p.u. specified active injection
Qinjection = (Qg(:) - Qd(:)) / Sbase; % p.u. specified reactive injection

V     = V0(:);           % initialize voltage magnitudes (p.u.)
delta = zeros(n,1);      % initialize voltage angles (rad), slack reference at 0

% convenience index sets for Jacobian construction and state updates
angle_ix = setdiff((1:n).', is);   % indices for unknown angles (non-slack buses)
vmag_ix  = ipq;                    % indices for unknown magnitudes (PQ buses)

% calculate G and B parts of Y
G = real(Y);
B = imag(Y);

N = 0; % iteration counter
tstart = tic; % start timer

% Newton Raphson loop
for k = 1:maxiter
    % calculate current operating point powers
    Vcomplex = V .* exp(1j*delta);        % complex bus voltages
    S  = Vcomplex .* conj(Y * Vcomplex);        % calculated complex bus power injections
    P  = real(S);                   % calculated active power
    Q  = imag(S);                   % calculated reactive power

    % calculate power mismatches (Specified - Calculated)
    dP = Pinjection - P;
    dQ = Qinjection - Q;

    % create mismatch vector for relevant equations
    rhs = [dP(angle_ix); dQ(vmag_ix)]; % P mismatches for non-slack, Q mismatches for PQ

    % convergence
    if norm(rhs, inf) <= toler
        N = k - 1; % record last successful iteration number
        break      % exit loop if converged
    end

    % Jacobian matrix
    % calculate elements using current V and delta
    th  = delta - delta.';          % matrix of angle differences (delta_i - delta_j)
    costh = cos(th);
    sinth = sin(th);
    ViVj = V .* V.';                % matrix of voltage magnitude products (V_i * V_j)

    % Off-diagonal Jacobian block elements
    Hij =  ViVj .* ( G .* sinth - B .* costh );  % dP/dδ (i≠j)
    Nij =   V   .* ( G .* costh + B .* sinth );  % V*dP/dV (i≠j) - Note: scaled N block element
    Mij = -ViVj .* ( G .* costh + B .* sinth );  % dQ/dδ (i≠j)
    Lij =   V   .* ( G .* sinth - B .* costh );  % V*dQ/dV (i≠j) - Note: scaled L block element

    % assemble blocks and overwrite diagonal elements
    H = Hij; Nblk = Nij; M = Mij; L = Lij; % Start with off-diagonals
    for i = 1:n % Calculate and place diagonal elements
        H(i,i)    = -Q(i) - B(i,i)*V(i)^2;                 % Hii = dPi/dδi
        Nblk(i,i) =  P(i)/max(V(i),eps) + G(i,i)*V(i);    % Nii = Vi*dPi/dVi
        M(i,i)    =  P(i) - G(i,i)*V(i)^2;                 % Mii = dQi/dδi
        L(i,i)    =  Q(i)/max(V(i),eps) - B(i,i)*V(i);    % Lii = Vi*dQi/dVi
    end

    % extract submatrices for the reduced system Ax=b
    Hs = H(angle_ix, angle_ix);     % relevant H block
    Ns = Nblk(angle_ix, vmag_ix);   % relevant N block
    Ms = M(vmag_ix,  angle_ix);     % relevant M block
    Ls = L(vmag_ix,  vmag_ix);      % relevant L block

    J  = [Hs, Ns; Ms, Ls];      % assemble the reduced Jacobian

    % solve linearly
    dx = J \ rhs;               % solve J * dx = rhs for dx using efficient backslash operator

    % update
    nang = numel(angle_ix); % number of unknown angles
    delta(angle_ix) = delta(angle_ix) + dx(1:nang);             % update angles
    V(vmag_ix)      = V(vmag_ix)      + dx(nang+1:end);         % update PQ voltage magnitudes

    N = k; % update iteration counter
end % end iteration loop

time = toc(tstart); % stop timer

% recalculate powers with final V and delta
Vcomplex    = V .* exp(1j*delta);
Scalc = Vcomplex .* conj(Y * Vcomplex);
Pcalc = real(Scalc);
Qcalc = imag(Scalc);

% calculate slack bus power and PV bus reactive power
Pslack = Pcalc(is) * Sbase + Pd(is);     % Slack active power in MW
Qgv = Qcalc(ipv) * Sbase + Qd(ipv);   % PV reactive power in Mvar

% ensure slack angle is exactly zero (reference)
delta = delta - delta(is);
delta(is) = 0;

end