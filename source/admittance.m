function Y = admittance(nfrom, nto, r, x, b)
% Deniz Temurcu 261089503
% This function computes the admittance matrix of a power system given line
% parameters

% Our inputs:
% nfrom is the vector of sending-end bus indices
% nto is the vector of receiving-end bus indices
% r is the vector of line resistances 
% x is the vector of line reactances 
% b is the vector of line susceptances 

% Our output:
% Y is the bus admittance matrix of size nbus (square matrix)

% --- Input size check ---
if ~( length(nfrom)==length(nto) && length(nto)==length(r) && ...
      length(r)==length(x) && length(x)==length(b) )
    error('Input vectors nfrom, nto, r, x, and b must all be the same length.');
end


nbus = max(max(nfrom), max(nto));% get the number of buses
len = length(r);
Y = zeros(nbus, nbus); % initialize as 0s matrix

for k = 1:len
        y = 1/(r(k) + 1i*x(k));      % calculate series admittance
        Y(nfrom(k),nto(k)) = Y(nfrom(k),nto(k)) - y; % subtract admittance from off-diagonal terms
        Y(nto(k),nfrom(k)) = Y(nto(k),nfrom(k)) - y; % should be symmetrical
        Y(nfrom(k),nfrom(k)) = Y(nfrom(k),nfrom(k)) + y + 1i*b(k)/2; % add admittances and half of shunt susceptances for diagonal terms
        Y(nto(k),nto(k))   = Y(nto(k),nto(k))   + y + 1i*b(k)/2;  
end
end
