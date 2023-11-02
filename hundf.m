function T = hundf(Nmax)
% Numbers of distinct electron arrangements of electrons in an orbital that satisfy
% the Pauli exclusion principle 
% INPUT:
%                                               s, p, d,  f
% Nmax   - maximum number of orbital electrons, 2, 6, 10, 14,... 
%          (4*n+2 for natural n including 0)
% OUTPUT:
% T      - numbers of distinct electron arrangements
% Based on
% https://pubs.aip.org/aip/adv/article/13/10/105308/2915332/The-second-law-of-infodynamics-and-its
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 17.10.2023 first working version, certainly not optimal
% v2: 02.11.2023 with sequence found by Andrzej Tomski

if nargin < 1
  error 'Wrong number of arguments in hundf.';
end

if Nmax<0
  error 'Nmax must be positive integer.';
end

for N=1:Nmax/2
    T(N) = floor((N+2)/2) * ( floor((N+2)/2) + 1 )/2;
end
for N=Nmax/2:Nmax
    T(N) = floor((Nmax-N+2)/2) * ( floor((Nmax-N+2)/2) + 1 )/2;
end

return    
