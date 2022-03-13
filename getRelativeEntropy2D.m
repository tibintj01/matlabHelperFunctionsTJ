function [D] = getRelativeEntropy2D(p,q)
    D=NaN;
    if(abs(sum(p(:))-1)>0.001)
        disp('p is not a probability distribution')
        return
    end
   
    
     if(abs(sum(q(:))-1)>0.001)
        disp('q is not a probability distribution')
        return
     end
    
    D=0;
    
    for i=1:size(p,1)
        for j=1:size(p,2)
            if(p(i,j)~=0 && q(i,j)~=0 )    
                D=D+p(i,j)*log(p(i,j)/q(i,j))/log(2);
            end
        end
    end
    
    if(D==0)
        D=NaN;
    end


