function D=fill_consensusclustering(A,gamma,reps,thresh,method,mode)

D=A;
[M,~]=size(D);
D(1:M+1:end)=0;

% fill in missing values with randomly selected existing values

D_nc=fill_missingvals(D,reps);

% community detection of ensemble of networks

switch method
    
    case 'Louvain'

    BF=zeros(M,M,reps);

    parfor ncidx=1:reps,

        Z=D_nc{ncidx};   
    
        V  = [1:size(Z,1)];
        Q0 = -1; Q1 = 0;            
        while Q1-Q0>1e-5;           
            Q0 = Q1;                
            [V, Q1] = community_louvain(Z,gamma,V);
        end
      
        BF(:,:,ncidx)=clustersol_representation(V);
      
    end

    case 'Infomap'
        
    BF=zeros(M,M,reps);

    for ncidx=1:reps,

        Z=D_nc{ncidx};   
        V=infomap(Z);
        
        BF(:,:,ncidx)=clustersol_representation(V);
      
    end             
        
end

% mean matrix

D=mean(BF,3);

if strcmp(mode,'binary') && strcmp(method,'Louvain'),
       
    D=double(D>=thresh);
    
elseif strcmp(mode,'binary') && strcmp(method,'Infomap'),
    
    D=D.*double(D>=thresh);
    
end
       
    %  switch method
    %         
    %     case 'Louvain'
    %     
    %      D=double(D>=thresh);
    %         
    %     case 'Infomap'
    %             
    %      D=D.*double(D>=thresh);
    %     
    %  end
    %     
    %  elseif strcmp(mode,'weighted'),
    %                
    %      D=D.*double(D>=thresh);
    % 
    %  end

% setting diagonal elements to zero

D(1:M+1:end)=0;

end