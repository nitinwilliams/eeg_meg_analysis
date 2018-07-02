function [PS]=partsim(PART1,PART2)

% Reference
% Ben-Hur A, Elisseeff A, and Guyon I. A stability based method for discovering structure in
% clustered data. Pacific Symposium on Biocomputing, 7:617, 2002.

% computing similarity between matrices
    
B1=clustersol_representation(PART1);  
B2=clustersol_representation(PART2);  

% dot-product
B12=B1.*B2;
    
% similarity terms
    
L12=sum(B12(:));
L11=sum(B1(:));
L22=sum(B2(:));

% partition similarity
    
PS=L12/sqrt(L11*L22);

end