clear all

% Matlab script to calculate and draw chemical elements entropy and spin multiplicity following Hund and Aufbau rule.
% Based on
% https://pubs.aip.org/aip/adv/article/13/10/105308/2915332/The-second-law-of-infodynamics-and-its
% https://www.preprints.org/manuscript/202310.1112

% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 16.10.2023 1st working version
% v2: 17.10.2023 code simplification
% v3: 18.10.2023 code simplification, spin multiplicities
% v3: 18.10.2023 automatic value generation (TODO: Aufbau exceptions)
% v4: 28.10.2023 A167268 and A216607 sequences

nmax=19; %19 = regular table  2900; 
n=1:nmax;

%otab = 4*mod(-n, round(sqrt(n))) + 2; % OEIS A167268 sequence
otab = 4*(floor( ceil( sqrt(4*n) ).^2/4 ) - n) + 2; % OEIS A216607 sequence

Z = 1; % Zmin
for k=1:length(otab)
    for N = 1:otab(k)
        Normax = otab(k);
        EL(Z, 1) = Z;
        EL(Z, 3) = (k-1)*log2(2);      % core entropy
        if N <= Normax/2
            EL(Z, 2) = N/2;             % spin multiplicity
        else
            EL(Z, 2) = (Normax-N)/2; % spin multiplicity
            EL(Z, 3) = EL(Z, 3) + log2(N) - (Normax/2)*log2(Normax/2)/N - (N-Normax/2)*log2(N-Normax/2)/N;  % surplus entropy
        end            
        EL(Z, 4) = otab(k);            % entropy multiplicity
        Z = Z + 1;
    end
end

s2 = log2(2);
d7 = log2(7) - 5*log2(5)/7 - (7-5)*log2(7-5)/7;  % surplus entropy
d8 = log2(8) - 5*log2(5)/8 - (8-5)*log2(8-5)/8;  % surplus entropy
d9 = log2(9) - 5*log2(5)/9 - (9-5)*log2(9-5)/9;  % surplus entropy

ELexc=[
%Z,	'S_act''conf_act'	'S_auf;			'conf_auf'      'cmt'      
24,		6, 5*s2;%			4; % h spin		6*s2,	    'exc',		
28,		2, 5*s2+d9;%		2;				6*s2+d8,	'exc',		
29,		1, 6*s2;%			1;				6*s2+d9,	'exc',		
41,		5, 8*s2;%			3; % h spin		9*s2,	    'exc',		
42,		6, 8*s2;%			4; % h spin		9*s2,	    'exc',		
44,		4, 8*s2+d7;%		4;				9*s2+d6,	'exc',		
45,		3, 8*s2+d8;%		3;				9*s2+d7,	'exc',		
46,		0, 9*s2;%			2; % l spin		9*s2+d8,	'exc',		
47,		1, 9*s2;%			1;				9*s2+d9,	'exc',		
57,		1, 12*s2;%			1;				12*s2,	    'exc',		
58,		2, 12*s2;%			2;				12*s2,	    'exc',		
64,		8, 12*s2;%			6; % h spin		12*s2+f8,	'exc',		
78,		2, 12*s2+d9;%		2;				13*s2+d8,	'exc',		
79,		1, 13*s2;%			1;				13*s2+d9,	'exc',		
89,		1, 16*s2;%			1;				16*s2,	    'exc',		
90,		2, 16*s2;%			2;				16*s2,	    'exc',		
91,		3, 16*s2;%			3;				16*s2,	    'exc',		
92,		4, 16*s2;%			4;				16*s2,	    'exc',		
93,		5, 16*s2;%			5;				16*s2,	    'exc',		
96,		8, 16*s2;%			6; % h spin		16*s2+f8,	'exc',		
103,	1, 17*s2];%			1;				17*s2,	    'exc',		

EL2=EL;

idx=1;
for k=1:size(EL, 1)
    if EL(k,1) == ELexc(idx,1) 
        EL2(k,2) = ELexc(idx,2)/2; %spin
        EL2(k,3) = ELexc(idx,3); %entropy        
        idx = idx+1;
    end
    if idx > size(ELexc, 1)
        break
    end
end

drawspin    = 1;
drawentropy = 1;

if drawspin
    sf = figure
    hold on
    grid on
    linew=1;
    plot(EL2(:, 1), EL2(:, 2), 'g', 'LineWidth', linew) %exceptions
    plot(EL(:, 1), EL(:, 2), 'r', 'LineWidth', linew)
    set(gca,'FontName', 'Times New Roman')
    set(gca,'FontSize', 12)
    xlabel('Atomic number')
    ylabel('Spin multiplicity')
    axis([2 max(EL2(:,1)) 0  max(EL2(:,2))])    

    % draw noble gases lines
    ngline = 0;
    for k=1:length(otab)
        if otab(k) == 2
            line([ngline ngline], [0 max(EL(:,2))], 'Color',[0 0 0], 'LineStyle', '-.');    
        end
        ngline = ngline + otab(k);            
    end

    % draw exceptions
     for k=1:size(ELexc,1)
         plot(ELexc(k, 1), ELexc(k, 2)/2, 'gx','LineWidth',linew)        
     end
    
    % resize figure
    %rect = [left, bottom, width, height]
    rect = get(sf, 'OuterPosition')
    rect(4) = rect(4)*.6;
    set(sf, 'OuterPosition', rect)
end
if drawentropy
    ef = figure
    hold on
    grid on
    linew=1;
    plot(EL2(:, 1),EL2(:, 3), 'g','LineWidth',linew) %exceptions    
    plot(EL(:, 1),EL(:, 3), 'r','LineWidth',linew)    
    set(gca,'FontName', 'Times New Roman')
    set(gca,'FontSize', 12)
    xlabel('Atomic number')
    ylabel('Shannon entropy (bits)')
    axis([1 max(EL(:,1)) 0  max(EL(:,3))])

    % draw noble gases lines
    ngline = 0;
    for k=1:length(otab)
        if otab(k) == 2
            line([ngline ngline], [0 max(EL(:,3))], 'Color',[0 0 0], 'LineStyle', '-.');
        end
        ngline = ngline + otab(k);        
    end

    % draw exceptions
     for k=1:size(ELexc,1)
         plot(ELexc(k, 1), ELexc(k, 3), 'gx','LineWidth',linew)        
     end
    % draw approximations
     %for k=1.5:0.01:1.6   
     %for k=1.542:0.001:1.545        
        k=1.6;
        plot(EL(:, 1), (EL(:, 1)-1).^(1/k), 'g','LineWidth',linew)
        ff = (EL(:, 1)-1).^log2(2);
        plot(EL(:, 1), ff, 'b','LineWidth',linew)        
     %end
end
