clear all

% Matlab script to calculate periodic table of elements
% Based on
% https://pubs.aip.org/aip/adv/article/13/10/105308/2915332/The-second-law-of-infodynamics-and-its
% https://www.preprints.org/manuscript/202310.1112

% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 22.10.2023 1st version

trows=3; % 3 => regular table, 80 => Z= 721762 used in the paper
cols = trows+2;
for r=1:trows
    for c=1:cols
        if c==1 
            ntab(r, c) = 2;
            continue
        end
        if c==2
            ntab(r, c) = 4*r+2-1;
            continue
        end
        ntab(r, c) = 0;
    end        
end
%ntab

mtab = zeros(trows-1, cols);
for r=1:trows-1
    for c=1:cols-r-1
        mtab(r,c) = 4*c+2;    
    end
    mtab(r,cols-r-1) =1;
end
%mtab

mtabl = zeros(1, cols);
mtabl(1) = 1;
mtab = [mtab; mtabl];    
mtab = rot90(mtab, 2)
ntab = ntab+mtab;

% add hydrogen and helium
HHe = zeros(1,cols);
HHe(1)            = 1; 
HHe(size(ntab,2)) = 1; 
ntab = [HHe; ntab]

%tab = ntab(1:rows-1, :)
tab = ntab
%return

rows = size(tab, 1);
cols = size(tab, 2);

% 1. calculate the length of the bottom row
brl = sum(tab(rows,:))

rw = 1; % row width, vertical % DO NOT CHANGE

pointsH  = [0 0;     1   0; 0     rw; 1   rw]; % hydrogen
pointsHe = [brl-1 0; brl 0; brl-1 rw; brl rw]; % helium 

points  = [pointsH; pointsHe];
pointsb = [1 0]; % boundary points  

figure
hold on
grid on
labsize   = 6;
labexcolor=[1 0 0];
set(gca,'FontName', 'Times New Roman')
set(gca,'FontSize', 16) % for periodic table
%set(gca,'FontSize', 12) % otherwise
xlabel('Group')
ylabel('Period')

line([0 1], [0 0]); % H 
line([0 1], [1 1]);
line([0 0], [0 1]);
line([1 1], [0 1]);
text(0.5, 0.5, '1','FontName', 'Times New Roman', 'FontSize',labsize)
line([brl-1 brl  ], [0 0]); % He 
line([brl-1 brl  ], [1 1]);
line([brl-1 brl-1], [0 1]);
line([brl   brl  ], [0 1]);
text((2*brl-1)/2, 0.5, '2','FontName', 'Times New Roman', 'FontSize',labsize)

Zexc = [24 28 29 41 42 44 45 46 47 57 58 64 78 79 89 90 91 92 93 96 103];

Zlabel = 3;
for r=2:rows % no fckn hydrogen
    for run=1:2
        rl = 0; % row length, horizontal
        fzeros = 0; % found zeros
        ftime  = 1; % last sections flag
        for c=1:cols   
            if tab(r, c) ~= 0
                if ~fzeros % first sections => draw points, add zeros
                    y1 = rw;
                    y2 = rw+1;
                    x1 = rl;
                    x2 = rl+tab(r, c);
                    prl= rl;
                    rl = rl+tab(r, c);
                    lgt=rl-prl; 
                    ty = (2*rw-1)/2 + 1;
                    tx = (2*prl+1)/2;
                    for k=0:1:lgt-1
                        if ismember(Zlabel, Zexc)
                            text(tx, ty, num2str(Zlabel),'FontName', 'Times New Roman', 'FontSize',labsize, 'Color', labexcolor)
                        else
                            text(tx, ty, num2str(Zlabel),'FontName', 'Times New Roman', 'FontSize',labsize)
                        end
                        Zlabel = Zlabel+1;
                        tx = tx+1;
                    end
                    
                    points = [points; x1, y1];
                    points = [points; x2, y1];
                    points = [points; x1, y2];
                    points = [points; x2, y2];

                    line([x1 x2], [y1 y1]);                    
                    line([x1 x2], [y2 y2]);                                        
                    line([x1 x1], [y1 y2]);                                        
                    line([x2 x2], [y1 y2]);                                                            
                    
                    if tab(r, c) ~= tab(rows, c) 
                        rl = rl+tab(rows, c)-tab(r, c);
                    end
                else % last sections => add zeros, draw points 
                    if tab(r, c) ~= tab(rows, c) 
                        rl = rl+tab(rows, c)-tab(r, c);
                    end
                    y1 = rw;
                    y2 = rw+1;
                    x1 = rl;
                    x2 = rl+tab(r, c);
                    prl= rl;
                    rl = rl+tab(r, c);
                    lgt=rl-prl;                     
                    ty = (2*rw-1)/2 + 1;
                    tx = (2*prl+1)/2;
                    for k=0:1:lgt-1
                        if ismember(Zlabel, Zexc)
                            text(tx, ty, num2str(Zlabel),'FontName', 'Times New Roman', 'FontSize',labsize, 'Color', labexcolor)
                        else
                            text(tx, ty, num2str(Zlabel),'FontName', 'Times New Roman', 'FontSize',labsize)
                        end
                        Zlabel = Zlabel+1;
                        tx = tx+1;
                    end
                    
                    %if ftime % are we 1st time here?
                    if ftime && run==1 % are we 1st time here and this is the 1st run?
                        if r ~= 2
                            pointsb = [pointsb; x2-1, y2-1];                    
                        end                            
                        ftime = 0;
                    end    
                    points = [points; x1, y1];
                    points = [points; x2, y1];
                    points = [points; x1, y2];
                    points = [points; x2, y2];

                    line([x1 x2], [y1 y1]);                    
                    line([x1 x2], [y2 y2]);                                        
                    line([x1 x1], [y1 y2]);                                        
                    line([x2 x2], [y1 y2]);                                                            
                end                
            else
                %if ~fzeros % this is the first time here
                if ~fzeros && run==1% this is the first time here                    
                           % so the previous points are boundary points
                        %if r ~= 2
                           pointsb = [pointsb; x2, y2-1];
                        %end
                end
                fzeros = 1; % found zeros
                rl = rl+tab(rows, c);
            end
        end
        rw = rw+1;    
    end
end    
% The number of elements is Zmax = Zlabel-1
Zlabel-1

pointsb = [pointsb; 2*(2*trows+2)-1, 2*trows-1];
%pointsb = [pointsb; 2*(2*trows+2),   2*trows-1];
pointsb = [pointsb; brl-1, 0];

points = sortrows(points);

%axis([0 32 0  7])    % regular periodic table
axis([0 brl 0  rw])

%return

%for k=1:size(points, 1) 
%    plot(points(k,1), points(k,2), 'bo');
%end    

pointsb = sortrows(pointsb);

for k=1:size(pointsb, 1) 
    plot(pointsb(k,1), pointsb(k,2), 'go');
end

%plot(pointsb(:,1), pointsb(:,2), 'g');

return

%{
% generate periodic table multplicities
% this is the initial tab
% without the H He row
cols=rows+1;
itab = zeros(rows,cols); % 
for r=1:rows
    for c=1:r+1
        itab(r,c) = 4*(c-1)+2;
    end
end

% WOW!
% This code automatically adds hydrogen!
% That means that our dividing assumption is OK
for r=1:rows
    for c=1:cols
        if c==r
            itab(r,c)   = itab(r,c)-1;
            itab(r,c+1) = 1;
        end
    end        
end        

% This code takes initial tab
% and makes a known periodic table out of it
% without the H He row
tabtoswap = tab(1:rows, 2:cols);
for r=1:rows
    tabtoswap(r,:) = fliplr(tabtoswap(r,:));
end
tab(1:rows, 2:cols) = tabtoswap;
return
%}
