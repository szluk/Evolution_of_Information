clear all

% Matlab script to populate electron orbitals following Hund rule.
% Based on the second law of information dynamics disclosed in
% https://pubs.aip.org/aip/adv/article/13/10/105308/2915332/The-second-law-of-infodynamics-and-its
% It clearly demonstartes the evolution of information that takes place since Big-Bang. 

% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 16.10.2023

%                                                  s, p, d,  f
Nmax=6;    % maximum number of orbital electrons, 2, 6, 10, 14,... 
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
disp('create all possible populations')
for k=1:idx
  confN{k}
end  

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

block

idx1=0;
for k=1:idx
  if ~ismember(k, block)
    idx1 = idx1 + 1;
    confN1{idx1} = confN{k};
  end
end  

disp('after removing unphysical situations')
for k=1:idx1
  confN1{k}
end  
size(confN1)

% remove negative s
% move them from confN1 to confN2 using idx
idx = 0;
for k=1:idx1
  if sum(sum(confN1{k})) >= 0
    idx = idx + 1;
    confN2{idx} = confN1{k};
  end
end

disp('after removing negative s')
for k=1:idx
  confN2{k}
end  
size(confN2)

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
sz=size(HH)

% this is Melvin configuration
HH = sortrows(HH,3)
nel = 1;
cidx=0;
cile=0;    
k=1;
while k <= size(HH,1)
  if ismember(nel, HH(k,3) )
    cidx = cidx+1;
    k = k+1;
  else
    cile(nel) = cidx;
    nel = nel+1;
    cidx= 0;
  end
end
cile(length(cile)+1)=1;    
    
cile    
sum(cile)
    
% remove duplicate entries having the same spin and entropy
HH = unique(HH,'rows');
HH = sortrows(HH,3)

figure
hold on
grid on
nel = 1;
tidx=0;
k=1;
while k <= size(HH,1)
  if ismember(nel, HH(k,3) )
    tidx= tidx+1;  
    Htp(tidx, 1) = HH(k,1);  
    Htp(tidx, 2) = HH(k,2);
    k = k+1;
  else
    plot( Htp(:,1), Htp(:,2), 'g-' )
    plot( Htp(:,1), Htp(:,2), 'ro' )    
    ile(nel) = tidx;
    nel = nel+1;
    clear Htp;
    tidx= 0;
  end
end
plot( Htp(:,1), Htp(:,2), 'g-' ) %plot last set ????
plot( Htp(:,1), Htp(:,2), 'ro' ) %plot last set ????
set(gca,'FontName', 'Times New Roman')
set(gca,'FontSize', 12)
xlabel('Spin multiplicity')
ylabel('Shannon entropy (bits)')
%axis([0 1.5 0.5 log(2)]) % Nmax=6,nats
%axis([0 2.5 0.4 log(2)]) % Nmax=10,nats

%log2(Nk) - Nmax*(log2(Nmax/2))./(2*Nk) - (2*Nk-Nmax).*log2((2*Nk-Nmax)/2)./(2*Nk)
%log2(Nk) - Nmax*(log2(Nmax)-log2(2))./(2*Nk) - (2*Nk-Nmax).*log2((2*Nk-Nmax)/2)./(2*Nk)
%log2(Nk) - Nmax*log2(Nmax)./(2*Nk) - (2*Nk-Nmax).*log2(2*Nk-Nmax)./(2*Nk) + log2(2)
%log2(Nk) - Nmax*log2(Nmax)./(2*Nk) - log2(2*Nk-Nmax) + Nmax*log2(2*Nk-Nmax)./(2*Nk) + log2(2)

ile

return

N1s = [1;
       1];
N1p = [1 0 0;
       1 2 1;
       2 3 2;
       2 3 2;
       3 0 0;
       3 0 0];
N1d = [1 0 0 0 0 0;
       1 2 1 0 0 0;
       2 3 2 0 0 0;
       2 3 2 4 3 2;
       3 4 3 5 4 3;
       3 4 3 5 4 3;       
       4 5 4 0 0 0;
       4 5 4 0 0 0;       
       5 0 0 0 0 0;
       5 0 0 0 0 0];
N1f = [1 0 0 0 0 0 0 0 0 0;
       1 2 1 0 0 0 0 0 0 0;
       2 3 2 0 0 0 0 0 0 0;
       2 3 2 4 3 2 0 0 0 0;
       3 4 3 5 4 3 0 0 0 0;
       3 4 3 5 4 3 6 4 3 0;      
       4 5 4 6 5 4 7 6 5 4;
       4 5 4 6 5 4 7 6 5 4;
       5 6 5 7 6 5 0 0 0 0;
       5 6 5 7 6 5 0 0 0 0;
       6 7 6 0 0 0 0 0 0 0;
       6 7 6 0 0 0 0 0 0 0;       
       6 0 0 0 0 0 0 0 0 0;              
       7 0 0 0 0 0 0 0 0 0];

N1act = N1d;
