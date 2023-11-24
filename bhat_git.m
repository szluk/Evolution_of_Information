clear all

% black hole assembly index
% Based on
% https://novapublishers.com/shop/chapter-15-black-hole-horizons-as-patternless-binary-messages-and-markers-of-dimensionality/
% https://www.researchgate.net/publication/375884073_Assembly_Theory_of_Patternless_Binary_Messages_How_to_Assemble_a_Black_Hole
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 22.11.2023

Nmax=10000;
Nmax=30;
Nmax=500;
for NN=1:Nmax
    alb(NN) = bhat_alb(NN);
    aub(NN) = bhat_aub(NN);    
end    

figure
hold on
grid on

% Marshal stuff
%axis([0 Nmax 0 max(aub)])
%axis([0 max(aub)+15 0 max(aub)])
%line([0 max(aub)], [-1   max(aub)-1], 'Color',[0 0 0], 'LineStyle', '-.');

% -1 stuff
axis([0 20 -1 13])
line([0 25], [0   0], 'Color',[0 0 0]);
plot([0 1], [-1 0], 'g:')
set(gca,'XTick', 0:1:25)
set(gca,'YTick',-1:1:15)
set(gca,'YTick',[-1 0 1 3 5 7 9 11 13])

plot(1:Nmax, alb, 'r')
plot(1:Nmax, aub, 'g')

plot(1:Nmax, log2(1:Nmax), 'r-.')
line([6  max(aub)], [5  max(aub)-1], 'Color',[0 1 0], 'LineStyle', '-.');

set(gca,'FontName', 'Times New Roman')
set(gca,'FontSize', 12)
xlabel('{\itN}')
ylabel('{\ita_N}')

return

% Saggittarius alb, aub
Nmax=7.24E90
%Nmax=13
k1 = (-3 + sqrt(9+4*Nmax))/2
k2 = (-3 - sqrt(9+4*Nmax))/2

kc = ceil(k1)

alb = bhat_alb(Nmax)
aub = 4*k1 - 1

aub = 4*kc - 1
return
