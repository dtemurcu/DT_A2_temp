function [V, delta, Pslack, Qgv, N, time] = decpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
% Deniz Temurcu 261089503
% This function performs Decoupled Newton power flow

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
% toler  is the convergence tolerance on mismatches (p.u.)
% maxiter is the maximum number of iterations (scalar)

% Our outputs:
% V      is the vector of bus voltage magnitudes (p.u., n x 1)
% delta  is the vector of bus voltage angles (rad, n x 1) with delta(is)=0
% Psl    is the slack-bus active generation (MW)
% Qgv    is the vector of PV-bus reactive generations (Mvar), same order as ipv
% N      is the number of iterations performed (scalar)
% time   is the CPU time (s) of the DECPF loop (scalar)
%

% check sizes
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

V     = V0(:);          % initialize voltage magnitudes (p.u.)
delta = zeros(n,1);     % initialize voltage angles (rad)

% index sets for Jacobian construction and state updates
angle_ix = setdiff((1:n).', is);   % indices for unknown angles (non-slack buses)
vmag_ix  = ipq;                    % indices for unknown magnitudes (PQ buses)

% calculate G and B parts of Y
G = real(Y);
B = imag(Y);

N = 0; % iteration counter
tstart = tic; % start timer

% Decoupled loop
for k = 1:maxiter
    % calculate current operating point powers
    Vcomplex = V .* exp(1j*delta);    % complex bus voltages
    S  = Vcomplex .* conj(Y * Vcomplex);    % calculated complex bus power injections
    P  = real(S);               % calculated active power
    Q  = imag(S);               % calculated reactive power

    % calculate power mismatches
    dP = Pinjection - P;
    dQ = Qinjection - Q;

    % convergence
    % Check ∞-norm of mismatches for the relevant subsystems
    if max([norm(dP(angle_ix),inf), norm(dQ(vmag_ix),inf)]) <= toler
        N = k - 1; % record last successful iteration number
        break      % exit loop if converged
    end

    % Decoupled Jacobian
    th  = delta - delta.';      % matrix of angle differences (delta_i - delta_j)
    costh = cos(th);
    sinth = sin(th);
    ViVj = V .* V.';            % matrix of voltage magnitude products (V_i * V_j)

    % off-diagonal elements
    Hij =  ViVj .* ( G .* sinth - B .* costh ); % dP/dδ (i≠j)
    Lij =   V   .* ( G .* sinth - B .* costh ); % V*dQ/dV (i≠j) - Scaled L element

    % assemble H and L blocks, overwrite diagonals
    H = Hij;  Lblk = Lij;           % start from off-diagonals
    for i = 1:n
        H(i,i)    = -Q(i) - B(i,i)*V(i)^2;                 % Hii = dPi/dδi
        Lblk(i,i) =  Q(i)/max(V(i),eps) - B(i,i)*V(i);    % Lii = Vi*dQi/dVi
    end

    % extract submatrices for the reduced decoupled systems
    Hs = H(angle_ix, angle_ix);         % Relevant H block (angles)
    Ls = Lblk(vmag_ix, vmag_ix);        % Relevant L block (magnitudes at PQ buses)

    % solve linear system
    ddelta = Hs \ dP(angle_ix);      % Solve Hs * ddelta = dP
    dV     = Ls \ dQ(vmag_ix);      % Solve Ls * dV = dQ (I renamed the variable for clarity)

    % update
    delta(angle_ix) = delta(angle_ix) + ddelta;       % update angles
    V(vmag_ix)      = V(vmag_ix) + dV;                % update mag

    N = k; % update iteration counter
end % end iteration loop

time = toc(tstart); % stop timer

% recalculate powers with final V and delta
Vcomplex    = V .* exp(1j*delta);
Scalc = Vcomplex .* conj(Y * Vcomplex);
Pcalc = real(Scalc);
Qcalc = imag(Scalc);

% calculate slack bus power and PV bus reactive power
Pslack = Pcalc(is) * Sbase + Pd(is);   % slack active power in MW
Qgv = Qcalc(ipv) * Sbase + Qd(ipv); % PV reactive power in Mvar

% output
% wnsure slack angle is exactly zero (reference)
delta = delta - delta(is);
delta(is) = 0;


end
