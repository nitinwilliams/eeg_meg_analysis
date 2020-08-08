function [PLV_mag,PLV_ang]=phasecoupling_estimation(D,PC,mode1,mode2)

M=size(D,2);
PLV_mag=zeros(M);
PLV_ang=zeros(M);

switch mode1 

    case 'undirected'
    
        if strcmp(mode2,'bivariate')
        
            for ridx=1:M

                for cidx=ridx:M
           
                    Y=D(:,ridx);
                    X=D(:,cidx);
                    B=X\Y;
                    
                    PLV_mag(ridx,cidx)=abs(B);
                    PLV_ang(ridx,cidx)=angle(B);
                    
                end
                
            end           
              
            PLV_mag=PLV_mag+triu(PLV_mag,1)';
            PLV_ang=PLV_ang+triu(PLV_ang,1)';
            
        elseif strcmp(mode2,'multivariate')
        
            for ridx=1:M          
                    
                    indices=setdiff([1:M],ridx);
                    
                    Y=D(:,ridx);
                    X=D(:,indices);
                    B=X\Y;
                    
                    PLV_mag(ridx,indices)=abs(B);
                    PLV_ang(ridx,indices)=angle(B);
                                   
            end   
                
        end
    
    case 'directed'
        
        if strcmp(mode2,'bivariate')
            
            for ridx=1:M

                for cidx=ridx+1:M
           
                    Y=[D(1:end-1,ridx),D(1:end-1,cidx)];
                    X=[D(2:end,ridx),D(2:end,cidx)];
                    
                    % B=Y\X;
                    B=X\Y;
                    
                    PLV_mag(ridx,cidx)=abs(B(1,2));
                    PLV_ang(ridx,cidx)=angle(B(1,2));
                    
                    PLV_mag(cidx,ridx)=abs(B(2,1));
                    PLV_ang(cidx,ridx)=angle(B(2,1));
                    
                end
                
            end     
            
        elseif strcmp(mode2,'multivariate')
            
                    Y=D(1:end-1,:);
                    X=D(2:end,:);
                    
                    % B=Y\X;
                    B=X\Y;
                    B=PC*B*pinv(PC);
                    
                    PLV_mag=abs(B);
                    PLV_ang=angle(B);         
                  
        end
       
end

end