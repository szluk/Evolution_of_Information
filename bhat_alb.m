function alb=bhat_alb(NBH)
% Matlab function to generate lower bound on a BH assembly index
% (OEIS sequence A014701)
% INPUT:
% NBH    - length of the string (BH information capacity)
% OUTPUT:
% alb    - lower bound on a BH assembly index
% Based on
% https://novapublishers.com/shop/chapter-15-black-hole-horizons-as-patternless-binary-messages-and-markers-of-dimensionality/
% https://www.researchgate.net/publication/375884073_Assembly_Theory_of_Patternless_Binary_Messages_How_to_Assemble_a_Black_Hole
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 18.11.2023 1st working version

alb = 0;
%%{
while NBH>1
   alb = alb + 1;
   if floor(NBH/2) ==  NBH/2 % even NBH
        NBH=NBH/2;
   else
        NBH=NBH-1;           
   end
end    
%%}

%{
first_step = true;           % first step flag
while NBH > 0
    for k=1:NBH
        if 2^k > NBH         % 2^k busted NBH
            break
        end
    end
    NBH = NBH - 2^(k-1);
    if first_step
        alb = alb+k-1;            
        first_step = false;
    else
        alb = alb+1;        
    end
end
%}