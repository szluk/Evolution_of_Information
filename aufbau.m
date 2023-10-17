clear all

% Matlab script to calculate and draw chemical elements' entropy following Hund and Aufbau rule.
% Based on
% https://pubs.aip.org/aip/adv/article/13/10/105308/2915332/The-second-law-of-infodynamics-and-its
% https://www.researchgate.net/publication/374752685_Shannon_Entropy_of_Chemical_Elements

% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 16.10.2023 1st working version
% v2: 17.10.2023 code simplification

%-1 in 'auf_conf' means that 'auf_conf' is the same as 'act_conf'
EL={ 
%nm, Z,   H_act_conf,                                                 H_auf_conf, comments,     act_conf#TODO#, auf_conf#TODO#
'H' ,1 ,  0,                                                          -1,      'hydrogen_has_zero_entropy';
'He',2 ,  log2(2),                                                    -1,      '1st_break';
'Li',3 ,  log2(2),                                                    -1,      'NONE';
'Be',4 ,  2*log2(2),                                                  -1,      'NONE';
'B', 5 ,  2*log2(2),                                                  -1,      'NONE';
'C', 6 ,  2*log2(2),                                                  -1,      'NONE';
'N', 7 ,  2*log2(2),                                                  -1,      'NONE';
'O', 8 ,  2*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4),   -1,      'NONE';
'F', 9 ,  2*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5),   -1,      'NONE';
'Ne',10,  3*log2(2),                                                  -1,      'NONE';
'Na',11,  3*log2(2),                                                  -1,      'NONE';
'Mg',12,  4*log2(2),                                                  -1,      'NONE';
'Al',13,  4*log2(2),                                                  -1,      'NONE';
'Si',14,  4*log2(2),                                                  -1,      'NONE';
'P', 15,  4*log2(2),                                                  -1,      'NONE';
'Si',16,  4*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4),   -1,      'NONE';
'Cl',17,  4*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5),   -1,      'NONE';
'Ar',18,  5*log2(2),                                                  -1,      '2nd_break';
'K', 19,  5*log2(2),                                                  -1,      'NONE';
'Ca',20,  6*log2(2),                                                  -1,      'NONE';
'Sc',21,  6*log2(2),                                                  -1,      'NONE';
'Ti',22,  6*log2(2),                                                  -1,      'NONE';
'V', 23,  6*log2(2),                                                  -1,      'NONE';
'Cr',24,  5*log2(2),                                                  6*log2(2), 'Cr_exception';
'Mn',25,  6*log2(2),                                                  -1,      'NONE';
'Fe',26,  6*log2(2)+(log2(6) - (5/6)*log2(5) - (6-5)*log2(6-5)/6),    -1,      'NONE';
'Co',27,  6*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7),    -1,      'NONE';
'Ni',28, 5*log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9),   6*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8), 'Ni disputed exception'; % [Ar] 3d8 4s2 % [Ar] 3d9 4s1 Nickel atomic
'Cu',29,  6*log2(2),                                                  6*log2(2)+(log2(9) - (5/9)*log2(5) - (9-5)*log2(9-5)/9),  'Cu_exception';% [Ar]+3d10+4s1
'Zn',30,  7*log2(2),                                                  -1,      'NONE';
'Ga',31,  7*log2(2),                                                  -1,      'NONE';
'Ge',32,  7*log2(2),                                                  -1,      'NONE';
'As',33,  7*log2(2),                                                  -1,      'NONE';
'Se',34,  7*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4),   -1,      'NONE';
'Br',35,  7*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5),   -1,      'NONE';
'Kr',36,  8*log2(2),                                                  -1,      '3rd_break';
'Rb',37,  8*log2(2),                                                  -1,      'NONE'; 
'Sr',38,  9*log2(2),                                                  -1,      'NONE';
'Y', 39,  9*log2(2),                                                  -1,      'NONE';
'Zr',40,  9*log2(2),                                                  -1,      'NONE';
'Nb',41,  9*log2(2),                                                  8*log2(2)+log2(2),      'Nb exception';
'Mo',42,  9*log2(2),                                                  8*log2(2)+log2(2),      'Mo exception';
'Tc',43,  9*log2(2),                                                  -1,      'NONE';
'Ru',44,  8*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7),    8*log2(2)+log2(2)+(log2(6) - (5/6) *log2(5) - (6-5) *log2(6-5)/6),      'Ru exception';
'Rh',45,  8*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8),    8*log2(2)+log2(2)+(log2(7) - (5/7) *log2(5) - (7-5) *log2(7-5)/7),      'Rh exception';
'Pd',46,  9*log2(2),                                                  8*log2(2)+log2(2)+(log2(8) - (5/8) *log2(5) - (8-5) *log2(8-5)/8),      'Pd exception';
'Ag',47,  9*log2(2),                                                  8*log2(2)+log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9),      'Ag exception';
'Cd',48,  10*log2(2),   -1,      'NONE';
'In',49,  10*log2(2),   -1,      'NONE';
'Sn',50,  10*log2(2),   -1,      'NONE';
'Sb',51,  10*log2(2),   -1,      'NONE';
'Te',52,  10*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4),                                                -1,      'NONE';
'I', 53,  10*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5),                                                -1,      'NONE';
'Xe',54,  11*log2(2),   -1,      '4th break';
'Cs',55,  11*log2(2),   -1,      'NONE';
'Ba',56,  12*log2(2),   -1,      'NONE';
'La',57,  12*log2(2),   11*log2(2)+log2(2),      'La exception';
'Ce',58,  12*log2(2),   11*log2(2)+log2(2),      'Ce exception';
'Pr',59,  12*log2(2),   -1,      'NONE';
'Nd',60,  12*log2(2),   -1,      'NONE';
'Pm',61,  12*log2(2),   -1,      'NONE';
'Sm',62,  12*log2(2),   -1,      'NONE';
'Eu',63,  12*log2(2),   -1,      'NONE';
'Gd',64,  12*log2(2),   11*log2(2)+log2(2)+(log2(8) - (7/8) *log2(7) - (8-7) *log2(8-7)/8),      'Gd exception';
'Tb',65,  12*log2(2)+(log2(9) - (7/9) *log2(7) - (9-7) *log2(9-7)/9),   -1,      'NONE';
'Dy',66,  12*log2(2)+(log2(10)- (7/10)*log2(7) - (10-7)*log2(10-7)/10),   -1,      'NONE';
'Ho',67,  12*log2(2)+(log2(11)- (7/11)*log2(7) - (11-7)*log2(11-7)/11),   -1,      'NONE';
'Er',68,  12*log2(2)+(log2(12)- (7/12)*log2(7) - (12-7)*log2(12-7)/12),   -1,      'NONE';
'Tm',69,  12*log2(2)+(log2(13)- (7/13)*log2(7) - (13-7)*log2(13-7)/13),   -1,      'NONE';
'Yb',70,  13*log2(2),   -1,      'NONE';
'Lu',71,  13*log2(2),   -1,      'NONE';
'Hf',72,  13*log2(2),   -1,      'NONE';
'Ta',73,  13*log2(2),   -1,      'NONE';
'W', 74,  13*log2(2),   -1,      'NONE';
'Re',75,  13*log2(2),   -1,      'NONE';
'Os',76,  13*log2(2)+(log2(6) - (5/6)*log2(5) - (6-5)*log2(6-5)/6),   -1,      'NONE';
'Ir',77,  13*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7),   -1,      'NONE';
'Pt',78,  13*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8),   11*log2(2)+log2(2)+log2(2)+(log2(8) - (5/8) *log2(5) - (8-5) *log2(8-5)/8),      'Pt exception';
'Au',79,  13*log2(2),   11*log2(2)+log2(2)+log2(2)+(log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9),      'Au exception';
'Hg',80,  14*log2(2),   -1,      'NONE';
'Tl',81,  14*log2(2),   -1,      'NONE';
'Pb',82,  14*log2(2),   -1,      'NONE';
'Bi',83,  14*log2(2),   -1,      'NONE';
'Po',84,  14*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4),   -1,      'NONE';
'At',85,  14*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5),   -1,      'NONE';
'Rn',86,  15*log2(2),   -1,      '5th break';
'Fr',87,  15*log2(2),   -1,      'NONE';
'Ra',88,  16*log2(2),   -1,      'NONE';
'Ac',89,  16*log2(2),   15*log2(2)+log2(2),      'Ac exception';
'Th',90,  16*log2(2),   15*log2(2)+log2(2),      'Th exception';
'Pa',91,  16*log2(2),   15*log2(2)+log2(2),      'Pa exception';
'U', 92   16*log2(2),   15*log2(2)+log2(2),      'U exception';
'Np',93,  16*log2(2),   15*log2(2)+log2(2),      'Np exception';
'Pu',94,  16*log2(2),   -1,      'NONE';
'Am',95,  16*log2(2),   -1,      'NONE';
'Cm',96,  16*log2(2),   15*log2(2)+log2(2)+(log2(9) - (7/9) *log2(7) - (9-7) *log2(9-7)/9),      'Cm exception';
'Bk',97,  16*log2(2)+(log2(9) - (7/9) *log2(7) - (9-7) *log2(9-7)/9),   -1,      'NONE';
'Cf',98,  16*log2(2)+(log2(10)- (7/10)*log2(7) - (10-7)*log2(10-7)/10),   -1,      'NONE';
'Es',99,  16*log2(2)+(log2(11)- (7/11)*log2(7) - (11-7)*log2(11-7)/11),   -1,      'NONE';
'Fm',100, 16*log2(2)+(log2(12)- (7/12)*log2(7) - (12-7)*log2(12-7)/12),   -1,      'NONE';
'Md',101, 16*log2(2)+(log2(13)- (7/13)*log2(7) - (13-7)*log2(13-7)/13),   -1,      'NONE';
'No',102, 17*log2(2),   -1,      'NONE';
'Lf',103, 17*log2(2),   15*log2(2)+log2(2)+log2(2),      'Lf exception';
'Rf',104, 17*log2(2),   -1,      'NONE';
'Db',105, 17*log2(2),   -1,      'NONE';
'Sg',106, 17*log2(2),   -1,      'NONE';
'Bh',107, 17*log2(2),   -1,      'NONE';
'Hs',108, 17*log2(2)+(log2(6) - (5/6)*log2(5) - (6-5)*log2(6-5)/6),   -1,      'NONE';
'Mt',109, 17*log2(2)+(log2(7) - (5/7)*log2(5) - (7-5)*log2(7-5)/7),   -1,      'NONE';
'Ds',110, 17*log2(2)+(log2(8) - (5/8)*log2(5) - (8-5)*log2(8-5)/8),   17*log2(2) + (log2(9) - (5/9) *log2(5) - (9-5) *log2(9-5)/9),      'NONE';
'Rg',111, 17*log2(2)+(log2(9) - (5/9)*log2(5) - (9-5)*log2(9-5)/9),   17*log2(2) + log2(2),      'NONE';
'Cn',112, 18*log2(2),   -1,      'NONE';
'Nh',113, 18*log2(2),   -1,      'NONE';
'Fl',114, 18*log2(2),   -1,      'NONE';
'Mc',115, 18*log2(2),   -1,      'NONE';
'Lv',116, 18*log2(2)+(log2(4) - (3/4)*log2(3) - (4-3).*log2(4-3)/4),   -1,      'NONE';
'Ts',117, 18*log2(2)+(log2(5) - (3/5)*log2(3) - (5-3).*log2(5-3)/5),   -1,      'NONE';
'Og',118, 19*log2(2),   -1,      'NONE'
};

for k=1:size(EL, 1)
  ZZ(k) = EL{k,2};
  HH(k) = EL{k,3};  
  if EL{k,4} == -1
    HHH(k)= EL{k,3};
    dfc(k)=0;
  else
    HHH(k)= EL{k,4};
    dfc(k)=1;    
  end    
end  

figure
hold on
grid on
linew=1;
plot(ZZ,HHH, 'r','LineWidth',linew)
plot(ZZ,HH, 'g','LineWidth',linew)
%plot(NiLineZ,NiLineen, 'g:','LineWidth',linew)
set(gca,'FontName', 'Times New Roman')
set(gca,'FontSize', 9)
xlabel('Atomic number')
ylabel('Shannon entropy (bits)')
axis([0 118 0 19]) % bits

for k=1:size(EL, 1)
  if HH(k) ~= HHH(k) % draw Aufbau rule violations with different entropy
     plot(ZZ(k), HHH(k),'rx','LineWidth',linew)
     plot(ZZ(k), HH(k),'gx','LineWidth',linew)
  else              % draw Aufbau rule violations with the same entropy
    if dfc(k) ~= 0
      plot(ZZ(k), HH(k),'rx','LineWidth',linew)
      plot(ZZ(k), HH(k),'b+','LineWidth',linew)  
    end  
  end
end  

return

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

% simplified formulas verification
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
