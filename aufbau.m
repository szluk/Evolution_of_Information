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
% v3: 18.10.2023 code simplification, looking for patterns...

s2 =   log2(2);
p4 = 2*log2(2)                 -log2(3^(3/4));
p5 =           -log2(2^(2/5))  -log2(3^(3/5))  +log2(5);
d6 =   log2(2)                                 -log2(5^(5/6)) +log2(3);
d7 =           -log2(2^(2/7))                  -log2(5^(5/7)) +log2(7);
d8 = 3*log2(2)                 -log2(3^(3/8))  -log2(5^(5/8)); 
d9 =           -log2(2^(8/9))                  -log2(5^(5/9)) +log2(9);
f8 = 3*log2(2)                                                          -log2(7^(7/8));
f9 =           -log2(2^(2/9))                                           -log2(7^(7/9))  +log2(9); 
f10=   log2(2)                 -log2(3^(3/10))                          -log2(7^(7/10)) +log2(5);
f11=           -log2(2^(8/11))                                          -log2(7^(7/11)) +log2(11);
f12= 2*log2(2)                                 -log2(5^(5/12))          -log2(7^(7/12)) +log2(3);
f13=           -log2(2^(6/13)) -log2(3^(6/13))                          -log2(7^(7/13)) +log2(13);

p4c = 2*log2(2) - (3/4)*log2(3);
p5c =   log2(5) - (3/5)*log2(3) - 2*log2(2)/5;

d6c = log2(6) -(5/6) *log2(5);
d7c = log2(7) -(5/7) *log2(5) - 2*log2(2)/7;
d8c = 3*log2(2) -(5/8) *log2(5) - 3*log2(3)/8;
d9c = log2(9) -(5/9) *log2(5) - 4*log2(4)/9;

f8c = 3*log2(2) -(7/8) *log2(7);
f9c = log2(9) -(7/9) *log2(7) -2*log2(2)/9;
f10c= log2(2) +log2(5)- (7/10)*log2(7) -3*log2(3)/10;
f11c= log2(11) -(7/11)*log2(7) -4*log2(4)/11;
f12c= log2(12) -(7/12)*log2(7) -5*log2(5)/12;
f13c= log2(13) -(7/13)*log2(7) -6*log2(6)/13;

EL={
%Z,	'conf_act'		'conf_auf';
1,	0,				0,	        'none';
2,	s2,			    s2,	        '1st break';
3,	s2,	  		    s2,	        'none';
4,	2*s2,		    2*s2,	    'none';
5,	2*s2,			2*s2,	    'none';
6,	2*s2,			2*s2,	    'none';
7,	2*s2,			2*s2,	    'none';
8,	2*s2+p4,		2*s2+p4,	'none';
9,	2*s2+p5,		2*s2+p5,	'none';
10,	3*s2,			3*s2,	    '2nd break';
11,	3*s2,			3*s2,	    'none';
12,	4*s2,			4*s2,	    'none';
13,	4*s2,			4*s2,	    'none';
14,	4*s2,			4*s2,	    'none';
15,	4*s2,			4*s2,	    'none';
16,	4*s2+p4,		4*s2+p4,	'none';
17,	4*s2+p5,		4*s2+p5,	'none';
18,	5*s2,			5*s2,	    '3rd break';
19,	5*s2,			5*s2,	    'none';
20,	6*s2,			6*s2,	    'none';
21,	6*s2,			6*s2,	    'none';
22,	6*s2,			6*s2,	    'none';
23,	6*s2,			6*s2,	    'none';
24,	5*s2,			6*s2,	    'exc'; %-s2
25,	6*s2,			6*s2,	    'none';
26,	6*s2+d6,		6*s2+d6,	'none';
27,	6*s2+d7,		6*s2+d7,	'none';
28,	5*s2+d9,		6*s2+d8,	'exc';%disp -d8+d9-s2
29,	6*s2,			6*s2+d9,	'exc';%-d9
30,	7*s2,			7*s2,	    'none';
31,	7*s2,			7*s2,	    'none';
32,	7*s2,			7*s2,	    'none';
33,	7*s2,			7*s2,	    'none';
34,	7*s2+p4,		7*s2+p4,	'none';
35,	7*s2+p5,		7*s2+p5,	'none';
36,	8*s2,			8*s2,	    '4th break';
37,	8*s2,			8*s2,	    'none';
38,	9*s2,			9*s2,	    'none';
39,	9*s2,			9*s2,	    'none';
40,	9*s2,			9*s2,	    'none';
41,	8*s2,			9*s2,	    'exc';%s2
42,	8*s2,			9*s2,	    'exc';%s2
43,	9*s2,			9*s2,	    'none';
44,	8*s2+d7,		9*s2+d6,	'exc';
45,	8*s2+d8,		9*s2+d7,	'exc';
46,	9*s2,			9*s2+d8,	'exc';
47,	9*s2,			9*s2+d9,	'exc';
48,	10*s2,			10*s2,	    'none';
49,	10*s2,			10*s2,	    'none';
50,	10*s2,			10*s2,	    'none';
51,	10*s2,			10*s2,	    'none';
52,	10*s2+p4,		10*s2+p4,	'none';
53,	10*s2+p5,		10*s2+p5,	'none';
54,	11*s2,			11*s2,	    '5th break';
55,	11*s2,			11*s2,	    'none';
56,	12*s2,			12*s2,	    'none';
57,	12*s2,			12*s2,	    'exc';
58,	12*s2,			12*s2,	    'exc';
59,	12*s2,			12*s2,	    'none';
60,	12*s2,			12*s2,	    'none';
61,	12*s2,			12*s2,	    'none';
62,	12*s2,			12*s2,	    'none';
63,	12*s2,			12*s2,	    'none';
64,	12*s2,			12*s2+f8,	'exc';
65,	12*s2+f9,		12*s2+f9,	'none';
66,	12*s2+f10,		12*s2+f10,	'none';
67,	12*s2+f11,		12*s2+f11,	'none';
68,	12*s2+f12,		12*s2+f12,	'none';
69,	12*s2+f13,		12*s2+f13,	'none';
70,	13*s2,			13*s2,	    'none';
71,	13*s2,			13*s2,	    'none';
72,	13*s2,			13*s2,	    'none';
73,	13*s2,			13*s2,	    'none';
74,	13*s2,			13*s2,	    'none';
75,	13*s2,			13*s2,	    'none';
76,	13*s2+d6,		13*s2+d6,	'none';
77,	13*s2+d7,		13*s2+d7,	'none';
78,	12*s2+d9,		13*s2+d8,	'exc';
79,	13*s2,			13*s2+d9,	'exc';
80,	14*s2,			14*s2,	    'none';
81,	14*s2,			14*s2,	    'none';
82,	14*s2,			14*s2,	    'none';
83,	14*s2,			14*s2,	    'none';
84,	14*s2+p4,		14*s2+p4,	'none';
85,	14*s2+p5,		14*s2+p5,	'none';
86,	15*s2,			15*s2,	    '6th break';
87,	15*s2,			15*s2,	    'none';
88,	16*s2,			16*s2,	    'none';
89,	16*s2,			16*s2,	    'exc';
90,	16*s2,			16*s2,	    'exc';
91,	16*s2,			16*s2,	    'exc';
92,	16*s2,			16*s2,	    'exc';
93,	16*s2,			16*s2,	    'exc';
94,	16*s2,			16*s2,	    'none';
95,	16*s2,			16*s2,	    'none';
96,	16*s2,			16*s2+f8,	'exc';
97,	16*s2+f9,		16*s2+f9,	'none';
98,	16*s2+f10,		16*s2+f10,	'none';
99,	16*s2+f11,		16*s2+f11,	'none';
100,16*s2+f12,		16*s2+f12,	'none';
101,16*s2+f13,		16*s2+f13,	'none';
102,17*s2,			17*s2,	    'none';
103,17*s2,			17*s2,	    'exc';
104,17*s2,			17*s2,	    'none';
105,17*s2,			17*s2,	    'predict';
106,17*s2,			17*s2,	    'predict';
107,17*s2,			17*s2,	    'predict';
108,17*s2+d6,		17*s2+d6,	'predict';
109,17*s2+d7,		17*s2+d7,	'predict';
110,17*s2+d8,		17*s2+d8,	'predict';
111,17*s2+d9,		17*s2+d9,	'predict';
112,18*s2,		    18*s2,	    'predict';
113,18*s2,		    18*s2,	    'predict';
114,18*s2,		    18*s2,	    'predict';
115,18*s2,		    18*s2,	    'predict';
116,18*s2+p4,		18*s2+p4,	'predict';
117,18*s2+p5,		18*s2+p5,	'predict';
118,19*s2,		    19*s2,	    'predict'};

for k=1:size(EL, 1)
  ZZ(k) = EL{k,1};
  HH(k) = EL{k,2};  
  HHH(k)= EL{k,3};
end  

figure
hold on
grid on
linew=1;
plot(ZZ,HHH, 'r','LineWidth',linew)
plot(ZZ,HH, 'g','LineWidth',linew)
%plot(NiLineZ,NiLineen, 'g:','LineWidth',linew)
set(gca,'FontName', 'Times New Roman')
set(gca,'FontSize', 6)
xlabel('Atomic number')
ylabel('Shannon entropy (bits)')
%axis([0 118 0 19*log2(2)]) % bits
axis([0 108 0 18*log2(2)]) % bits

for k=1:size(EL, 1)
    if strcmp(EL{k,4}, 'exc')
        plot(ZZ(k), HHH(k),'rx','LineWidth',linew)
        plot(ZZ(k), HH(k),'gx','LineWidth',linew)        
        if   HH(k) == HHH(k)
          plot(ZZ(k), HH(k),'b+','LineWidth',linew)  
        end    
    end
    if k == 46
        plot(ZZ(k), HHH(k),'ro','LineWidth',linew)
        plot(ZZ(k), HH(k), 'go','LineWidth',linew)        
    end    
end  

line([2   2], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');
line([10 10], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');
line([18 18], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');
line([36 36], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');
line([54 54], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');
line([86 86], [0 19], 'Color',[0 0 0], 'LineStyle', '-.');

line([0 118], [9 9], 'Color',[0 0 0], 'LineStyle', '-.');

% simplified formulas verification
% Z S_act S_Auf H_act H_Auf
TabEx= [24  6 4 5 *log2(2)                         6*log2(2);
        28  2 2  log2(2^(37/9) * 3^2 * 5^(-5/9))   log2(2^9 * 3^(-3/8) * 5^(-5/8) );
        29	1 1 6 *log2(2)                         log2(2^(46/9) * 3^2 * 5^(-5/9));
        41  5 3 8 *log2(2)                         9*log2(2);
        42  6 4 8 *log2(2)                         9*log2(2);
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
