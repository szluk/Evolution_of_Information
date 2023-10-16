clear all

% Matlab script to populate calculate chemical elements' entropy following Hund and Aufbau rule.
% Based on the second law of information dynamics disclosed in
% https://pubs.aip.org/aip/adv/article/13/10/105308/2915332/The-second-law-of-infodynamics-and-its
% https://www.researchgate.net/publication/374752685_Shannon_Entropy_of_Chemical_Elements
% It clearly demonstartes the evolution of information that takes place since Big-Bang. 

% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 16.10.2023

%s
N = 2;
Nk=1:N;

% s1=0
% log2(2)=log2(2)

%log2(Nk) - (N/2)*log2(N/2) - (Nk-N/2).*log2(Nk-N/2)./Nk
%Ns=2
%Np=6
%Nd=10
%Nf=14

EL={
%nm Z  conf
'H' ,1 ,  0;
'He',2 ,  log2(2);   % 1st break
'Li',3 ,  log2(2);
'Be',4 ,  2*log2(2);
'B', 5 ,  2*log2(2);
'C', 6 ,  2*log2(2);
'N', 7 ,  2*log2(2);
'O', 8 ,  2*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4);
'F', 9 ,  2*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5);
'Ne',10,  3*log2(2);
'Na',11,  3*log2(2);
'Mg',12,  4*log2(2);
'Al',13,  4*log2(2);
'Si',14,  4*log2(2);
'P', 15,  4*log2(2);
'Si',16,  4*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4);
'Cl',17,  4*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5);
'Ar',18,  5*log2(2); % 2nd break
'K', 19,  5*log2(2);
'Ca',20,  6*log2(2);
'Sc',21,  6*log2(2);
'Ti',22,  6*log2(2);
'V', 23,  6*log2(2);
%'Cr',24, 6*log2(2); % Chromium as should be
'Cr',24,  5*log2(2); % Chromium exception
'Mn',25,  6*log2(2);
'Fe',26,  6*log2(2)+(log2(6) - (5/6)*log2(5) - (6-5)*log2(6-5)/6);
'Co',27,  6*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7);
'Ni',28,  6*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8); % [Ar] 3d8 4s2 Nickel Aufbau
%'Ni',28, 5*log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9); % [Ar] 3d9 4s1 Nickel atomic
%'Cu',29, 6*log2(2)+(log2(9) - (5/9)*log2(5) - (9-5)*log2(9-5)/9);  % Copper as should be
'Cu',29,  6*log2(2); % Copper exception [Ar]+3d10+4s1
'Zn',30,  7*log2(2);
'Ga',31,  7*log2(2);
'Ge',32,  7*log2(2);
'As',33,  7*log2(2);
'Se',34,  7*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4);
'Br',35,  7*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5);
'Kr',36,  8*log2(2); % 3rd break
'Rb',37,  8*log2(2); 
'Sr',38,  9*log2(2);
'Y', 39,  9*log2(2);
'Zr',40,  9*log2(2);
'Nb',41,  9*log2(2); % exc
'Mo',42,  9*log2(2); % exc
'Tc',43,  9*log2(2);
'Ru',44,  8*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7); % exc
'Rh',45,  8*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8); % exc
'Pd',46,  9*log2(2); % exc
'Ag',47,  9*log2(2); % exc
'Cd',48,  10*log2(2);
'In',49,  10*log2(2);
'Sn',50,  10*log2(2);
'Sb',51,  10*log2(2);
'Te',52,  10*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4);
'I', 53,  10*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5);
'Xe',54,  11*log2(2); % 4th break
'Cs',55,  11*log2(2);
'Ba',56,  12*log2(2);
'La',57,  12*log2(2); % exc
'Ce',58,  12*log2(2); % exc
'Pr',59,  12*log2(2);
'Nd',60,  12*log2(2);
'Pm',61,  12*log2(2);
'Sm',62,  12*log2(2);
'Eu',63,  12*log2(2);
'Gd',64,  12*log2(2); % exc
'Tb',65,  12*log2(2)+(log2(9) - (7/9) *log2(7) - (9-7) *log2(9-7)/9);
'Dy',66,  12*log2(2)+(log2(10)- (7/10)*log2(7) - (10-7)*log2(10-7)/10);
'Ho',67,  12*log2(2)+(log2(11)- (7/11)*log2(7) - (11-7)*log2(11-7)/11);
'Er',68,  12*log2(2)+(log2(12)- (7/12)*log2(7) - (12-7)*log2(12-7)/12);
'Tm',69,  12*log2(2)+(log2(13)- (7/13)*log2(7) - (13-7)*log2(13-7)/13);
'Yb',70,  13*log2(2);
'Lu',71,  13*log2(2);
'Hf',72,  13*log2(2);
'Ta',73,  13*log2(2);
'W', 74,  13*log2(2);
'Re',75,  13*log2(2);
'Os',76,  13*log2(2)+(log2(6) - (5/6)*log2(5) - (6-5)*log2(6-5)/6);
'Ir',77,  13*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7);
'Pt',78,  13*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8); % exc
'Au',79,  13*log2(2); % exc
'Hg',80,  14*log2(2);
'Tl',81,  14*log2(2);
'Pb',82,  14*log2(2);
'Bi',83,  14*log2(2);
'Po',84,  14*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4);
'At',85,  14*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5);
'Rn',86,  15*log2(2); % 5th break
'Fr',87,  15*log2(2);
'Ra',88,  16*log2(2);
'Ac',89,  16*log2(2); % exc
'Th',90,  16*log2(2); % exc
'Pa',91,  16*log2(2); % exc
'U', 92   16*log2(2); % exc
'Np',93,  16*log2(2); % exc
'Pu',94,  16*log2(2);
'Am',95,  16*log2(2);
'Cm',96,  16*log2(2); % exc
'Bk',97,  16*log2(2)+(log2(9) - (7/9) *log2(7) - (9-7) *log2(9-7)/9);
'Cf',98,  16*log2(2)+(log2(10)- (7/10)*log2(7) - (10-7)*log2(10-7)/10);
'Es',99,  16*log2(2)+(log2(11)- (7/11)*log2(7) - (11-7)*log2(11-7)/11);
'Fm',100, 16*log2(2)+(log2(12)- (7/12)*log2(7) - (12-7)*log2(12-7)/12);
'Md',101, 16*log2(2)+(log2(13)- (7/13)*log2(7) - (13-7)*log2(13-7)/13);
'No',102, 17*log2(2);
'Lf',103, 17*log2(2); % exc
'Rf',104, 17*log2(2);
'Db',105, 17*log2(2);
'Sg',106, 17*log2(2);
'Bh',107, 17*log2(2);
'Hs',108, 17*log2(2)+(log2(6) - (5/6)*log2(5) - (6-5)*log2(6-5)/6);
'Mt',109, 17*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7);
'Ds',110, 17*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8);
'Rg',111, 17*log2(2)+(log2(9) - (5/9)*log2(5) - (9-5)*log2(9-5)/9);
'Cn',112, 18*log2(2);
'Nh',113, 18*log2(2);
'Fl',114, 18*log2(2);
'Mc',115, 18*log2(2);
'Lv',116, 18*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4);
'Ts',117, 18*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5);
'Og',118, 19*log2(2)};

for k=1:size(EL, 1)
  ZZ(k) = EL{k,2};
  HH(k) = EL{k,3};  
end  

%s1 = 0;                                                            Ss1=1
%s2 = log2(2);                                                      Ss2=0

%p1 to p3 = 0                                                       Spk=k     
%p4 = (log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4)                Sp4=2
%p5 = (log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5)                Sp5=1 
%p6 = (log2(6) - (3/6)*log2(3) - (6-3).*log2(6-3)/6 = log2(2)       Sp6=0

%d1 to d5 = 0                                                       Sdk=k
%d6 = (log2(6) - (5/6) *log2(5) - (6-5) *log2(6-5)/6)               Sd6=4
%d7 = (log2(7) - (5/7) *log2(5) - (7-5) *log2(7-5)/7)               Sd7=3  
%d8 = (log2(8) - (5/8) *log2(5) - (8-5) *log2(8-5)/8)               Sd8=2
%d9 = (log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9)               Sd9=1
%d10= (log2(10)- (5/10)*log2(5) - (10-5)*log2(10-5)/10) = log2(2)   Sd9=0

%f1 to f7 = 0                                                       Sfk=k 
%f8 = (log2(8) - (7/8) *log2(7) - (8-7) *log2(8-7)/8)               Sf8=6
%f9 = (log2(9) - (7/9) *log2(7) - (9-7) *log2(9-7)/9)               Sf8=5
%f10= (log2(10)- (7/10)*log2(7) - (10-7)*log2(10-7)/10)             Sf8=4 
%f11= (log2(11)- (7/11)*log2(7) - (11-7)*log2(11-7)/11)             Sf8=3
%f12= (log2(12)- (7/12)*log2(7) - (12-7)*log2(12-7)/12)             Sf8=2
%f13= (log2(13)- (7/13)*log2(7) - (13-7)*log2(13-7)/13)             Sf8=1 
%f14= (log2(14)- (7/14)*log2(7) - (14-7)*log2(14-7)/14) = log2(2)   Sf8=0

% Chromium exception 24
%[Ar]+3d3+4s2, [Ar]+3d5+4s1, [Ar]+3d5+4s2
%instead of
%[Ar]+3d3+4s2, [Ar]+3d4+4s2, [Ar]+3d5+4s2
% Cr_exc = [Ar]+3d4+4s2 = 5*log2(2) + 0 + log2(2) % Chromium as should be
Cr_exc = 6*log2(2);

%28 Nickel exception [Ar]+3d9+4s1
Ni_exc = 6*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8);

% Copper exception 29
%[Ar]+3d8+4s2, [Ar]+3d10+4s1, [Ar]+3d10+4s2
%instead of
%[Ar]+3d8+4s2, [Ar]+3d9+4s2, [Ar]+3d10+4s2
%Cu_exc = [Ar]+3d9+4s2 = 5*log2(2) + (log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9) + log2(2) % Copper as should be
Cu_exc = 6*log2(2)+(log2(9) - (5/9)*log2(5) - (9-5)*log2(9-5)/9);  

%41[Kr]+5s2+4d3
Nb_exc = 8*log2(2)+log2(2);
%42[Kr]+5s2+4d4
Mo_exc = 8*log2(2)+log2(2);
%44[Kr]+5s2+4d6
Ru_exc = 8*log2(2)+log2(2)+(log2(6) - (5/6) *log2(5) - (6-5) *log2(6-5)/6);
%45[Kr]+5s2+4d7
Rh_exc = 8*log2(2)+log2(2)+(log2(7) - (5/7) *log2(5) - (7-5) *log2(7-5)/7);
%46[Kr]+5s2+4d8
Pd_exc = 8*log2(2)+log2(2)+(log2(8) - (5/8) *log2(5) - (8-5) *log2(8-5)/8);
%47[Kr]+5s2+4d9
Ag_exc = 8*log2(2)+log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9);
%57[Xe]+6s2+4f1
La_exc = 11*log2(2)+log2(2)+0;
%58[Xe]+6s2+4f2
Ce_exc = 11*log2(2)+log2(2)+0;
%64[Xe]+6s2+4f8
Gd_exc = 11*log2(2)+log2(2)+(log2(8) - (7/8) *log2(7) - (8-7) *log2(8-7)/8);
%78[Xe]+6s2+4f14+5d8
Pt_exc = 11*log2(2)+log2(2)+log2(2)+(log2(8) - (5/8) *log2(5) - (8-5) *log2(8-5)/8);
%79[Xe]+6s2+4f14+5d9
Au_exc = 11*log2(2)+log2(2)+log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9);
%89[Rn]+7s2+5f1
Ac_exc = 15*log2(2)+log2(2);
%90[Rn]+7s2+5f2
Th_exc = 15*log2(2)+log2(2);
%91[Rn]+7s2+5f3
Pa_exc = 15*log2(2)+log2(2);
%92[Rn]+7s2+5f4
U_exc = 15*log2(2)+log2(2);
%93[Rn]+7s2+5f5
Np_exc = 15*log2(2)+log2(2);
%96[Rn]+7s2+5f8
Cm_exc = 15*log2(2)+log2(2)+(log2(9) - (7/9) *log2(7) - (9-7) *log2(9-7)/9);
%103[Rn]+7s2+5f14+6d1
Lf_exc = 15*log2(2)+log2(2)+log2(2);

%110
Ds_exc = 17*log2(2) + (log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9); % .	110	*[Rn]f14d9s1   
%111
Rg_exc = 17*log2(2) + log2(2); % 111	*[Rn]f14d10s1  

HHH = HH;
HHH(24)=Cr_exc;
%HHH(28)=Ni_exc;
HHH(29)=Cu_exc;
HHH(41)=Nb_exc;
HHH(42)=Mo_exc;
HHH(44)=Ru_exc;
HHH(45)=Rh_exc;
HHH(46)=Pd_exc;
HHH(47)=Ag_exc;
HHH(57)=La_exc;
HHH(58)=Ce_exc;
HHH(64)=Gd_exc;
HHH(78)=Pt_exc;
HHH(79)=Au_exc;
HHH(89)=Ac_exc;
HHH(90)=Th_exc;
HHH(91)=Pa_exc;
HHH(92)=U_exc;
HHH(93)=Np_exc;
HHH(96)=Cm_exc;
HHH(103)=Lf_exc;

% Nickel
%'Ni',28, 5*log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9); % [Ar] 3d9 4s1 Nickel atomic
Niat28 = 5*log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9);

NiLineZ = [27     28     29]; 
NiLineen= [HH(27) Niat28 Cu_exc];

figure
hold on
grid on
linew=1;
plot(ZZ,HHH, 'r','LineWidth',linew)
plot(ZZ,HH, 'g','LineWidth',linew)
plot(NiLineZ,NiLineen, 'g:','LineWidth',linew)
set(gca,'FontName', 'Times New Roman')
set(gca,'FontSize', 9)
xlabel('Atomic number')
ylabel('Shannon entropy (bits)')
axis([0 118 0 19]) % bits

%plot(110, Ds_exc, 'co')
%plot(111, Rg_exc, 'co')

%line([43   43], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');
%line([48   48], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');

plot(24, Cr_exc,'rx','LineWidth',linew)
plot(24, HH(24),'gx','LineWidth',linew)
%plot(28, Ni_exc,'rx','LineWidth',linew)
plot(29, Cu_exc,'rx','LineWidth',linew)
plot(29, HH(29),'gx','LineWidth',linew)
plot(41, Nb_exc,'rx','LineWidth',linew)
plot(42, Mo_exc,'rx','LineWidth',linew)
plot(44, Ru_exc,'rx','LineWidth',linew)
plot(44, HH(44),'gx','LineWidth',linew)
plot(45, Rh_exc,'rx','LineWidth',linew)
plot(45, HH(45),'gx','LineWidth',linew)
%plot(46, Pd_exc,'rx','LineWidth',linew)
plot(46, Pd_exc,'ro','LineWidth',linew)
plot(46, HH(46),'go','LineWidth',linew)
plot(47, Ag_exc,'rx','LineWidth',linew)
plot(47, HH(47),'gx','LineWidth',linew)
plot(57, La_exc,'rx','LineWidth',linew)
plot(58, Ce_exc,'rx','LineWidth',linew)
plot(64, Gd_exc,'rx','LineWidth',linew)
plot(64, HH(64),'gx','LineWidth',linew)
plot(78, Pt_exc,'rx','LineWidth',linew)
plot(79, Au_exc,'rx','LineWidth',linew)
plot(79, HH(79),'gx','LineWidth',linew)
plot(89, Ac_exc,'rx','LineWidth',linew)
plot(90, Th_exc,'rx','LineWidth',linew)
plot(91, Pa_exc,'rx','LineWidth',linew)
plot(92, U_exc,'rx','LineWidth' ,linew)
plot(93, Np_exc,'rx','LineWidth',linew)
plot(96, Cm_exc,'rx','LineWidth',linew)
plot(96, HH(96),'gx','LineWidth',linew)
plot(103,Lf_exc,'rx','LineWidth',linew)

%draw points with the same entropy
plot(41, Nb_exc,'b+','LineWidth',linew)
plot(42, Mo_exc,'b+','LineWidth',linew)
plot(57, La_exc,'b+','LineWidth',linew)
plot(58, Ce_exc,'b+','LineWidth',linew)
plot(78, Pt_exc,'b+','LineWidth',linew)
plot(89, Ac_exc,'b+','LineWidth',linew)
plot(90, Th_exc,'b+','LineWidth',linew)
plot(91, Pa_exc,'b+','LineWidth',linew)
plot(92, U_exc, 'b+','LineWidth',linew)
plot(93, Np_exc,'b+','LineWidth',linew)
plot(103,Lf_exc,'b+','LineWidth',linew)

% Z S_act S_Auf H_act H_Auf
TabEx= [24  6 4 5 *log2(2)                         6*log2(2);
        28  2 2  log2(2^(37/9) * 3^2 * 5^(-5/9))   log2(2^9 * 3^(-3/8) * 5^(-5/8) );
        29	1 1 6 *log2(2)                         log2(2^(46/9) * 3^2 * 5^(-5/9));
        41  5 3 9 *log2(2)                         9*log2(2);
        42  6 4 9 *log2(2)                         9*log2(2);
        44	4 4 log2(2^(54/7) * 5^(-5/7) * 7 )     log2(2^10 * 3 * 5^(-5/6));        
        45  3 3 log2(2^11 * 3^(-3/8) * 5^(-5/8))   log2(2^(61/7) * 5^(-5/7) * 7);        
        46  0 2 9 *log2(2)                         log2(2^12 * 3^(-3/8) * 5^(-5/8) );
        47  1 1 9 *log2(2)                         log2(2^(73/9) * 3^2 * 5^(-5/9));
        57  1 1 12*log2(2)                         12*log2(2);
        58  2 2 12*log2(2)                         12*log2(2);
        64	8 6 12*log2(2)                         log2(2^15 * 7^(-7/8));
        78  2 2 log2(2^16 * 3^(-3/8) *5^(-5/8))    log2(2^16 * 3^(-3/8) * 5^(-5/8));
        79	1 1 13*log2(2)                         log2(2^(109/9) *3^2 * 5^(-5/9));
        89  1 1 16*log2(2)                         16*log2(2);
        90  2 2 16*log2(2)                         16*log2(2);
        91  3 3 16*log2(2)                         16*log2(2);
        92  4 4 16*log2(2)                         16*log2(2);
        93  5 5 16*log2(2)                         16*log2(2);
        96	8 6 16*log2(2)                         log2(2^(142/9) * 3^2 * 7^(-7/9));
        103 1 1 17*log2(2)                         17*log2(2)];

%plot( TabEx(:,1),  TabEx(:,4), 'bo','LineWidth',linew) %act
%plot( TabEx(:,1),  TabEx(:,5), 'c+','LineWidth',linew) %Auf

%        28  2 2 (37/9)*log2(2)+2*log2(3) - (5/9) *log2(5)              9*log2(2) - (5/8)*log2(5) - (3/8)*log2(3);
%        29	1 1 6 *log2(2)                                             (46/9)*log2(2) + 2*log2(3) - (5/9)*log2(5);
%        44	4 4 (54/7)*log2(2) + log2(7) - (5/7)*log2(5)               10*log2(2) + log2(3) - (5/6)*log2(5);



