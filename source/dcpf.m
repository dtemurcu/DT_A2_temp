function [delta, Pslack, Pflow] = dcpf(nfrom, nto, x, is, Pg, Pd, Sbase)
% Deniz Temurcu 261089503
% This function performs DC (linear) power flow using branch reactances (series-only B')

% Our inputs:
% nfrom is the vector of sending-end bus indices (m x 1, 1-based)
% nto   is the vector of receiving-end bus indices (m x 1, 1-based)
% x     is the vector of branch series reactances (p.u. on Sbase, m x 1)
% is    is the slack bus index (scalar, 1-based)
% Pg    is the vector of active power generation (MW, n x 1)
% Pd    is the vector of active power demand (MW, n x 1)
% Sbase is the MVA base (scalar)

% Our outputs:
% delta is the vector of bus voltage angles (rad, n x 1) with delta(is)=0
% Psl   is the slack-bus active generation (MW)
% Pf    is the vector of branch active power flows From->To (MW, m x 1)
%

% size check
n = max(max(nfrom), max(nto)); % get the number of buses
m = length(nfrom); % get number of branches
if ~isscalar(is) || is < 1 || is > n
    error('Invalid slack index.');
end
if any([numel(Pg), numel(Pd)] ~= n)
    error('Pg and Pd must be length n.');
end
if any([length(nto), length(x)] ~= m)
    error('nfrom, nto, x must be the same length.');
end
nfrom = nfrom(:); nto = nto(:); x = x(:); Pg = Pg(:); Pd = Pd(:); % ensure column vectors

if any(x <= 0) % Reactance should be positive, zero causes division error
    warning('Input reactances `x` should be positive. Check input data.');
    x(x <= 0) = eps; % Replace non-positive with small epsilon to avoid errors
end

% B' is the susceptance matrix considering only series reactances (neglecting resistance and shunts)
w  = 1 ./ x; % vector of series susceptances (b = 1/x)

% b B' efficiently using sparse matrix operations
% off-diagonal elements: B'_ij = -b_ij
% diagonal elements: B'_ii = sum_{k connected to i} b_ik
Bp = sparse(nfrom, nto, -w, n, n);   % populate off-diagonal (-w) for i->j
Bp = Bp + sparse(nto, nfrom, -w, n, n);   % populate off-diagonal (-w) for j->i (symmetry)
diagonal_elements = -sum(Bp, 2);          % calculate diagonal elements as negative sum of off-diagonals in the row
Bp = Bp + spdiags(diagonal_elements, 0, n, n); % add diagonal elements

% solve for angles
angle_ix = setdiff((1:n).', is); % indices for unknown angles (non-slack buses)
rhs    = (Pg - Pd) / Sbase;    % p.u. net active power injections

delta  = zeros(n,1); % initialize voltage angles (rad)

% solve the reduced linear system: Bp_reduced * delta_unknown = rhs_known
Bp_reduced = Bp(angle_ix, angle_ix);
rhs_reduced = rhs(angle_ix);
delta(angle_ix) = Bp_reduced \ rhs_reduced; % solve using backslash operator

% outputs
% calculate all bus injections using the full Bp matrix and solved angles
Pinjection_calc = Bp * delta; % p.u. calculated injections

% slack bus power calculation
Pslack = Pinjection_calc(is) * Sbase + Pd(is); % Slack power in MW (Calculated Injection + Slack Demand)

% branch flow calculation
Pflow  = ((delta(nfrom) - delta(nto)) ./ x) * Sbase; % Branch flows in MW: (delta_i - delta_j) / x_ij * Sbase

% ensure slack angle is exactly zero (reference)
delta = delta - delta(is);
delta(is) = 0;

end