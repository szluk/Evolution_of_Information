function aub=bhat_aub(NBH)
% Matlab function to generate upper bound on a BH assembly index
% INPUT:
% NBH    - length of the string (BH information capacity)
% OUTPUT:
% aub    - upper bound on a BH assembly index
% Based on
% https://novapublishers.com/shop/chapter-15-black-hole-horizons-as-patternless-binary-messages-and-markers-of-dimensionality/
% https://www.researchgate.net/publication/375884073_Assembly_Theory_of_Patternless_Binary_Messages_How_to_Assemble_a_Black_Hole
% 
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 20.11.2023 1st working version

%N=0 a(0)=-1
% up, up, side 0, up, up, side 0
% up, up, side 1, up, up, side 1
% up, up, side 2, up, up, side 2
% up, up, side 3, up, up, side 3
% up, up, side 4, up, up, side 4

step=1;
run =1;
side=0;

NN = 0;
aub= -1;
while NN < NBH
    if step < 3
        NN = NN+1;
        aub= aub + 1;
    else % step==3
        for k=1:side
            if side > 0
                NN = NN+1;                    
            end
        end
        run = run+1;
        if run > 2
            run  = 1;
            side = side+1;
        end        
    end
    step = step+1;
    if step > 3
        step=1;
    end
end
