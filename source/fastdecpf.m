function [V, delta, Psl, Qgv, N, time] = fastdecpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
% Deniz Temurcu 261089503
% This function performs Fast–Decoupled Power Flow

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
% toler  is the convergence tolerance (p.u.)
% maxiter is the maximum number of FDPF iterations (scalar)

% Our outputs:
% V      is the vector of bus voltage magnitudes (p.u., n x 1)
% delta  is the vector of bus voltage angles (rad, n x 1) with delta(is)=0
% Psl    is the slack-bus active generation (MW)
% Qgv    is the vector of PV-bus reactive generations (Mvar), same order as ipv
% N      is the number of iterations performed (scalar)
% time   is the CPU time (s) of the FDPF loop (scalar)
%

% size check
n = size(Y,1); % get the number of buses
if size(Y,2) ~= n
    error('Y must be square.');
end
ipq = ipq(:); ipv = ipv(:); % ensure column vectors
if any([numel(Pg), numel(Qg), numel(Pd), numel(Qd), numel(V0)] ~= n)
    error('Pg,Qg,Pd,Qd,V0 must be length n.');
end
if any(ipq == is) || any(ipv == is)
    error('ipq and ipv must not include the slack bus.');
end
if ~isempty(intersect(ipq, ipv))
    error('ipq and ipv must be disjoint.');
end

% initlialize
Pinjection = (Pg(:) - Pd(:)) / Sbase; % p.u. specified active injection
Qinjection = (Qg(:) - Qd(:)) / Sbase; % p.u. specified reactive injection

V     = V0(:);          % initialize voltage magnitudes (p.u.)
delta = zeros(n,1);     % initialize voltage angles (rad), slack reference at 0

% index sets for Jacobian construction and state updates
angle_ix = setdiff((1:n).', is);   % indices for unknown angles (non-slack buses)
vmag_ix  = ipq;                    % indices for unknown magnitudes (PQ buses)

% build constant B' and B'' Matrices (BX method) ---
Bfull = imag(Y);                 % B = -imag(Y) as per standard definition, Bfull here is imag(Y)
Bseries = Bfull;                 % start with Bfull
Bdiag   = diag(Bseries);         % get diagonal elements (includes shunts)
Boff    = Bseries - diag(Bdiag); % get off-diagonal elements
% recalculate diagonals for B' using only series elements
% diagonal element B'_ii = -sum_{k!=i} B_ik (where B_ik = imag(Y_ik))
Bseries(1:n+1:end) = -sum(Boff, 2);

Bp  = -Bseries(angle_ix, angle_ix);         % B' matrix for angle subsystem (using series susceptances only)
Bpp = -Bfull(vmag_ix,  vmag_ix);            % B'' matrix for magnitude subsystem (using full susceptances)

% pre factorize using decomposition for efficiency
Fp = decomposition(Bp,  'lu');
Fq = decomposition(Bpp, 'lu');

N = 0; % iteration counter
tstart = tic; % start timer

% Fast Decoupled Loop
for k = 1:maxiter
    % calculate current operating point powers
    Vc = V .* exp(1j*delta);    % complex bus voltages
    S  = Vc .* conj(Y * Vc);    % calculated complex bus power injections
    P  = real(S);               % calculated active power
    Q  = imag(S);               % calculated reactive power

    % calculate power mismatches
    dP = Pinjection - P;
    dQ = Qinjection - Q;

    % convergence
    % check based on relevant mismatches
    if max([max(abs(dP(angle_ix))), max(abs(dQ(vmag_ix)))]) <= toler
        N = k - 1; % record last successful iteration number
        break      % exit loop if converged
    end

    % updates
    % Use voltage magnitude guard against division by zero if V approaches 0
    rhsP = dP(angle_ix) ./ max(V(angle_ix), 1e-6); % P mismatch divided by |V|
    rhsQ = dQ(vmag_ix)  ./ max(V(vmag_ix),  1e-6); % Q mismatch divided by |V|

    % solve linear system
    ddelta = Fp \ rhsP;          % Solve Bp * ddelta = rhsP
    dVpq   = Fq \ rhsQ;          % Solve Bpp * dV = rhsQ (Note: result is dV, not dV/V)

    % update
    delta(angle_ix) = delta(angle_ix) + ddelta; % update angles
    V(vmag_ix)      = V(vmag_ix)      + dVpq;   % update PQ voltage magnitudes

    N = k; % update iteration counter
end % end iteration loop

time = toc(tstart); % stop timer

% recalculate powers with final V and delta
Vc    = V .* exp(1j*delta);
Scalc = Vc .* conj(Y * Vc);
Pcalc = real(Scalc);
Qcalc = imag(Scalc);

% calculate slack bus power and PV bus reactive power
Psl = Pcalc(is) * Sbase + Pd(is);       % Slack active power in MW
Qgv = Qcalc(ipv) * Sbase + Qd(ipv);     % PV reactive power in Mvar

% output
% ensure slack angle is exactly zero (reference)
delta = delta - delta(is);
delta(is) = 0;

end
