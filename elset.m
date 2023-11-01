clear all

% Matlab script to calculate and draw chemical elements electron-set entropy following Hund and Aufbau rule.
% Based on
% https://www.preprints.org/manuscript/202310.1112
% https://www.researchgate.net/publication/374752685_Shannon_Entropy_of_Chemical_Elements

% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 01.11.2023 1st working version

nmax=19; %19 = regular table  2900; 
n=1:nmax;

otab1 = 4*(floor( ceil( sqrt(4*n) ).^2/4 ) - n) + 2;   % OEIS A167268 and A216607 sequences
otab2 = 4*(-floor(-(1/4)*ceil(sqrt(4*n)).^2) - n) + 2; % OEIS A167268 and A366932 sequences

Z = 1; % Zmin
Zmax = sum(otab1);
Nu= zeros(1, Zmax); % numbers of electron up
Nd= zeros(1, Zmax); % numbers of electron down
prevNu = 0;
prevNd = 0;    
for k=1:length(otab1)
    for N = 1:otab1(k)
        Normax = otab1(k);    
        if N <= Normax/2        
            Nu(Z) = prevNu + N;
        else
            Nu(Z) = prevNu + Normax/2;% + (Normax-N)+1;
        end
        Nd(Z) = Z-Nu(Z);
        if Z == 18
            Ar_Nu = Nu(Z);
        end
        if Z == 36
            Kr_Nu = Nu(Z);
        end
        if Z == 54
            Xe_Nu = Nu(Z);
        end
        if Z == 86
            Rn_Nu = Nu(Z);
        end
        Z = Z + 1;
    end
    prevNu = Nu(Z-1);
    prevNd = Nd(Z-1);    
end

H(1) = 0;
for Z=2:Zmax
   pu = Nu(Z)/Z; 
   pd = Nd(Z)/Z;
   H(Z) = pu*log2(1/pu) + pd*log2(1/pd);
end

%otab2 --------------------------

Z = 1; % Zmin
Zmax2 = sum(otab2);
Nu2= zeros(1, Zmax2); % numbers of electron up
Nd2= zeros(1, Zmax2); % numbers of electron down
prevNu = 0;
prevNd = 0;    
for k=1:length(otab2)
    for N = 1:otab2(k)
        Normax = otab2(k);    
        if N <= Normax/2        
            Nu2(Z) = prevNu + N;
        else
            Nu2(Z) = prevNu + Normax/2;% + (Normax-N)+1;
        end
        Nd2(Z) = Z-Nu2(Z);
        Z = Z + 1;
    end
    prevNu = Nu2(Z-1);
    prevNd = Nd2(Z-1);    
end

H2(1) = 0;
for Z=2:Zmax2
   pu = Nu2(Z)/Z; 
   pd = Nd2(Z)/Z;
   H2(Z) = pu*log2(1/pu) + pd*log2(1/pd);
end

%exceptions
Hex = H;
p24u = (Ar_Nu+6)/24; 	%Cr	[Ar]+3d5+4s1
p28u = (Ar_Nu+6)/28;    %Ni	[Ar]+3d9+4s1
p29u = (Ar_Nu+6)/29;    %Cu	[Ar]+3d10+4s1
p41u = (Kr_Nu+5)/41;    %Nb	[Kr]+4d4+5s1
p42u = (Kr_Nu+6)/42;	%Mo	[Kr]+4d5+5s1
p44u = (Kr_Nu+6)/44;	%Ru	[Kr]+4d7+5s1
p45u = (Kr_Nu+6)/45;	%Rh	[Kr]+4d8+5s1
p46u = (Kr_Nu+5)/46;	%Pd	[Kr]+4d10
p47u = (Kr_Nu+6)/47;	%Ag	[Kr]+4d10+5s1
p57u = (Xe_Nu+2)/57;	%La	[Xe]+5d1+6s2
p58u = (Xe_Nu+3)/58;	%Ce	[Xe]+4f1+5d1+6s2
p64u = (Xe_Nu+9)/64;	%Gd	[Xe]+4f7+5d1+6s2
p78u = (Xe_Nu+13)/78;	%Pt	[Xe]+4f14+5d9+6s1
p79u = (Xe_Nu+13)/79;	%Au	[Xe]+4f14+5d10+6s1
p89u = (Rn_Nu+2)/89;	%Ac	[Rn]+6d1+7s2
p90u = (Rn_Nu+3)/90;	%Th	[Rn]+6d2+7s2
p91u = (Rn_Nu+4)/91;	%Pa	[Rn]+5f2+6d1+7s2
p92u = (Rn_Nu+5)/92;	%U	[Rn]+5f3+6d1+7s2
p93u = (Rn_Nu+6)/93;	%Np	[Rn]+5f4+6d1+7s2
p96u = (Rn_Nu+9)/96;	%Cm	[Rn]+5f7+6d1+7s2
p103u= (Rn_Nu+9)/103;   %Lf	[Rn]+5f14+7s2+7p1

Hex(24) = p24u*log2(1/p24u) + (1-p24u)*log2(1/(1-p24u));
Hex(28) = p28u*log2(1/p28u) + (1-p28u)*log2(1/(1-p28u));
Hex(29) = p29u*log2(1/p29u) + (1-p29u)*log2(1/(1-p29u));
Hex(41) = p41u*log2(1/p41u) + (1-p41u)*log2(1/(1-p41u));
Hex(42) = p42u*log2(1/p42u) + (1-p42u)*log2(1/(1-p42u));
Hex(44) = p44u*log2(1/p44u) + (1-p44u)*log2(1/(1-p44u));
Hex(45) = p45u*log2(1/p45u) + (1-p45u)*log2(1/(1-p45u));
Hex(46) = p46u*log2(1/p46u) + (1-p46u)*log2(1/(1-p46u));
Hex(47) = p47u*log2(1/p47u) + (1-p47u)*log2(1/(1-p47u));
Hex(57) = p57u*log2(1/p57u) + (1-p57u)*log2(1/(1-p57u));
Hex(58) = p58u*log2(1/p58u) + (1-p58u)*log2(1/(1-p58u));
Hex(64) = p64u*log2(1/p64u) + (1-p64u)*log2(1/(1-p64u));
Hex(78) = p78u*log2(1/p78u) + (1-p78u)*log2(1/(1-p78u));
Hex(79) = p79u*log2(1/p79u) + (1-p79u)*log2(1/(1-p79u));
Hex(89) = p89u*log2(1/p89u) + (1-p89u)*log2(1/(1-p89u));
Hex(90) = p90u*log2(1/p90u) + (1-p90u)*log2(1/(1-p90u));
Hex(91) = p91u*log2(1/p91u) + (1-p91u)*log2(1/(1-p91u));
Hex(92) = p92u*log2(1/p92u) + (1-p92u)*log2(1/(1-p92u));
Hex(93) = p93u*log2(1/p93u) + (1-p93u)*log2(1/(1-p93u));
Hex(96) = p96u*log2(1/p96u) + (1-p96u)*log2(1/(1-p96u));
Hex(103)= p103u*log2(1/p103u) + (1-p103u)*log2(1/(1-p103u));

%    [24, 28, 29, 41, 42, 44, 45, 46, 47, 57, 58, 64, 78, 79, 89, 90, 91, 92, 93, 96, 103]
%     Cr, Ni, Cu, Nb, Mo, Ru, Rh, Pd, Ag, La, Ce, Gd, Pt, Au, Ac, Th, Pa, U,  Np, Cm, Lf, 
%$Z=\{24, 28, 29, 41, 42, 44, 45, 46, 47, 57, 58, 64, 78, 79, 89, 90, 91, 92, 93, 96, 103\}$ %Aufbau rule violators
Zex1=[    28, 29,         44, 45,     47, 57, 58,     78, 79, 89, 90, 91, 92, 93,     103];% violators entropies are the same
Zex2=[24,         41, 42,                         64,                             96     ];%violators spin multiplicity is higher
Zex3=[                            46                                                     ];%violator spin multiplicity is lower

sf = figure;
hold on
grid on
linew=1;
plot(1:Zmax, H,   'r', 'LineWidth', linew)
plot(1:Zmax, Hex, 'g', 'LineWidth', linew)
plot(1:Zmax2, H2, 'b', 'LineWidth', linew)
set(gca,'FontName', 'Times New Roman')
set(gca,'FontSize', 7)
xlabel('Atomic number')
ylabel('Shannon entropy (bits)')
                          
%for k=1:length(Zex1)
%	plot(Zex1(k), H(Zex1(k)),'rx','LineWidth',linew)   
%	plot(Zex1(k), Hex(Zex1(k)),'gx','LineWidth',linew)       
%end    
for k=1:length(Zex1)
	plot(Zex1(k), H(Zex1(k)),'b+','LineWidth',linew)   
	plot(Zex1(k), H(Zex1(k)),'rx','LineWidth',linew)   	
end    
for k=1:length(Zex2)
	plot(Zex2(k), H(Zex2(k)),'rx','LineWidth',linew)   
	plot(Zex2(k), Hex(Zex2(k)),'gx','LineWidth',linew)       
end    
plot(46, H(46),'ro','LineWidth',linew)
plot(46, Hex(46),'go','LineWidth',linew)

% draw noble gases lines
ngline = 0;
for k=1:length(otab1)
    if otab1(k) == 2
        line([ngline ngline], [0 1], 'Color',[0 0 0], 'LineStyle', '-.');
    end
    ngline = ngline + otab1(k);        
end

axis([1 118 0.979 1])

rect = get(sf, 'OuterPosition')
rect(4) = rect(4)*.6;
set(sf, 'OuterPosition', rect)

%{
"24	"		Cr	[Ar]+3d5+4s1
"28	"		Ni	[Ar]+3d9+4s1
"29	"		Cu	[Ar]+3d10+4s1
"41	"		Nb	[Kr]+4d4+5s1
"42	"		Mo	[Kr]+4d5+5s1
"44	"		Ru	[Kr]+4d7+5s1
"45	"		Rh	[Kr]+4d8+5s1
"46	"		Pd	[Kr]+4d10
"47	"		Ag	[Kr]+4d10+5s1
"57	"		La	[Xe]+5d1+6s2
"58	"		Ce	[Xe]+4f1+5d1+6s2
"64	"		Gd	[Xe]+4f7+5d1+6s2
"78	"		Pt	[Xe]+4f14+5d9+6s1
"79	"		Au	[Xe]+4f14+5d10+6s1
"89	"		Ac	[Rn]+6d1+7s2
"90	"		Th	[Rn]+6d2+7s2
"91	"		Pa	[Rn]+5f2+6d1+7s2
"92	"		U	[Rn]+5f3+6d1+7s2
"93	"		Np	[Rn]+5f4+6d1+7s2
"96	"		Cm	[Rn]+5f7+6d1+7s2
103		Lf	[Rn]+5f14+7s2+7p1
%}