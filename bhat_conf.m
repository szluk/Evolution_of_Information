function [snew s]=bhat_conf(NBH)
% Matlab function to generate balanced unique strigs of 0s and 1s
% (patternless black hole messages)
% INPUT:
% NBH    - length of the string (BH information capacity)
% OUTPUT:
% snew   - distinctive strings
% s      - binom(NBH, floor(NBH/2)) of strings with equal number of 0s and 1s
% Based on
% https://novapublishers.com/shop/chapter-15-black-hole-horizons-as-patternless-binary-messages-and-markers-of-dimensionality/
% https://www.researchgate.net/publication/375884073_Assembly_Theory_of_Patternless_Binary_Messages_How_to_Assemble_a_Black_Hole
% 
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 16.11.2023 1st working version
% v2: 18.11.2023 1st removed conf
% v3: 19.11.2023 1st new algo for reduction 

% 1. generate a configuration from all configurations and keep it if it satisfies
N1= floor(NBH/2);
idx = 1;
for k=1:2^NBH
    confk = dec2bin(k-1,NBH);
    for l=1:NBH    
        conf(l) = bin2dec(confk(1,l));
    end
    if sum( conf ) == N1
        s(idx,:) = conf;
        idx = idx+1;
    end
end    

% 2. reduce to distinctive ones = remove duplicates
r=0;
snew = s;
while 1
    clear ss
    clear matching_rows
    r = r+1; 
    if r > size(snew,1)    
        break
    end    
    ss = [snew snew];
    conf2find = snew(r,:);
    for k=1:size(snew,1)    
        foundr = strfind(ss(k,:), conf2find);
        if ~isempty(foundr)
            matching_rows(k) = 1;
        else
            matching_rows(k) = 0;
        end 
    end

    idx = 1;
    matching_row_inserted = false;
    clear ns    
    for k=1:size(snew,1)
        if matching_rows(k)
            if ~matching_row_inserted
                ns(idx,:) = snew(k,:);
                matching_row_inserted = true;
                idx = idx + 1;
            end
        else
                ns(idx,:) = snew(k,:);
                idx = idx + 1;                
        end
    end
    snew = ns;
end
