function dar = hundf(Nmax)
% Numbers of distinct electron arrangements of electrons in an orbital that satisfy
% the Pauli exclusion principle 
% INPUT:
%                                               s, p, d,  f
% Nmax   - maximum number of orbital electrons, 2, 6, 10, 14,... 
%          (4*n+2 for natural n including 0)
% OUTPUT:
% dar    - numbers of distinct electron arrangements
% Based on
% https://pubs.aip.org/aip/adv/article/13/10/105308/2915332/The-second-law-of-infodynamics-and-its
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 17.10.2023 first working version, certainly not optimal

if nargin < 1
  error 'Wrong number of arguments in hundf.';
end

if Nmax>18
  error 'Nmax>18 is not supported in this version of hundf - Calculation takes too long time.';
end

if Nmax<0
  error 'Nmax must be positive integer.';
end

for k=0:4
  nn = [2 6 10 14 18];
  if ~ismember(Nmax, nn)
     error 'Nmax must be in the form 4*n+2 for natural n including 0';
  end
end

nos = Nmax/2; % number of "steps"

% create blank configuration of a population left col (+1), right col (-1)
conf = zeros(nos, 2);
conf(:,1) = 1;
conf(:,2) =-1;

% create all binary sequences (strings)
lcols_str=dec2bin(2^nos-1:-1:0);

lcols = zeros(2^nos,nos);
rcols = lcols;

for r=1:2^nos
  for c=1:nos
    if str2num(lcols_str(r,c)) == 1
      lcols(r,c) = 1;
      rcols(r,c) =-1;      
    end
  end
end

%lcols
%rcols

% phase 0 - create all possible populations
idx=0;
for l=1:(2^nos)-1
  for r=1:2^nos
      idx = idx+1;
      confN{idx} = [lcols(l,:); rcols(r,:)];
  end
end  
%disp('create all possible populations')
%for k=1:idx
%  confN{k}
%end  

% remove unphysical situations
% move them from confN to confN1 using idx1
block=0; % for Nmax=2
for col = 1:nos-1;
  idx1=0;
  for k=1:idx
    % 1st check
    if ((confN{k}(1,col) == 0 && confN{k}(1,col+1) ==  1) ...
     || (confN{k}(2,col) == 0 && confN{k}(2,col+1) == -1))
        % 2nd check
        if ~((confN{k}(1,col) ==  1 && confN{k}(1,col+1) ==  0) ...
          && (confN{k}(2,col) ==  0 && confN{k}(2,col+1) == -1))
            idx1 = idx1 + 1;
            block(col,idx1) = k;
        end
    end
  end
end

block;

idx1=0;
for k=1:idx
  if ~ismember(k, block)
    idx1 = idx1 + 1;
    confN1{idx1} = confN{k};
  end
end  

%disp('after removing unphysical situations')
%for k=1:idx1
%  confN1{k}
%end  
%size(confN1)

% remove negative s
% move them from confN1 to confN2 using idx
idx = 0;
for k=1:idx1
  if sum(sum(confN1{k})) >= 0
    idx = idx + 1;
    confN2{idx} = confN1{k};
  end
end

%disp('after removing negative s')
%for k=1:idx
%  confN2{k}
%end  
%size(confN2)

% calculate Shanon entropy
for k=1:idx
  N =0;  % all electrons
  Np=0;  %  1 electrons
  Nn=0;  % -1 electrons  
  for s=1:nos
    if confN2{k}(1,s) == 1
      N = N+1;
      Np= Np+1;
    end   
    if confN2{k}(2,s) == -1
      N = N+1;
      Nn= Nn+1;
    end   
  end
  pp = Np/N;
  pn = Nn/N;  
  H  = 0;
  if pp ~= 0
    H = H + pp*log2(1/pp); %bits
    %H = H + pp*log(1/pp); %nats   
  end  
  if pn ~= 0   
    H = H + pn*log2(1/pn); %bits
    %H = H + pn*log(1/pn); %nats   
  end
  HH(k,1) = (Np-Nn)/2;% spin
  HH(k,2) = H;        % entropy
  HH(k,3) = N;        % number of electrons
end  
sz=size(HH);

HH = sortrows(HH,3);
nel = 1;
cidx=0;
dar=0;    
k=1;
while k <= size(HH,1)
  if ismember(nel, HH(k,3) )
    cidx = cidx+1;
    k = k+1;
  else
    dar(nel) = cidx;
    nel = nel+1;
    cidx= 0;
  end
end
dar(length(dar)+1)=1;    
    
return    
