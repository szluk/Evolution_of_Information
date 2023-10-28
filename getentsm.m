function [H S] = getentsm(Z)
% Matlab function to calculate entropy H and spin multiplicity S
% of an element having the atomic number Z
% following Hund and Aufbau rule.
% INPUT:
% Z      - atomic number
% OUTPUT:
% H      - Shannon entropy
% S      - spin multiplicity
% Based on
% https://www.preprints.org/manuscript/202310.1112
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 28.10.2023 1st working version

H = 0; % starting entropy
S = 0; % starting spin-multiplicity

kmax = ceil((Z-1)^log(2))+1; % bound for k
O = zeros(1,kmax);     % vector of orbitals {2, 6, 10, 14, 18,...}

% populate vector of orbitals
Zi= 0; % iterated Z
for k=1:kmax
    bk = floor( ceil( sqrt(4*k) ).^2/4 ) - k;
    ak = 4*bk + 2;       % nof orbital electrons 
    Oi = (ak-2)/4 + 1;   % index of the vector of orbitals to store    
    Zi = Zi + ak;    
    O(Oi) = O(Oi) + 1;            
    lastOi = Oi;         % last stored orbital
    Nleft = Zi - Z;      % nof not present electrons required to fill the orbital
    N     = ak - Nleft;  % nof electrons present at the orbital
    if Nleft >= 0
        break
    end
    H = H + log2(2);
end
O = nonzeros(O)';

if Nleft == 0; % this is pure-entropy-contributor, so
    H = H + log2(2); % add pure entropy 
else
    Normax = 4*(lastOi-1)+2;
    if N <= Normax/2
        S = N/2;            % spin multiplicity
    else
        S = (Normax-N)/2;   % spin multiplicity
        H = H + log2(N) - (Normax/2)*log2(Normax/2)/N - (N-Normax/2)*log2(N-Normax/2)/N;  % surplus entropy
    end
end    
return
