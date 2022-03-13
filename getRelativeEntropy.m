function [D] = getRelativeEntropy(p,q)
    D=0;
    if(abs(sum(p(:))-1)>0.001)
        disp('p is not a probability distribution')
        return
    end
   
     if(abs(sum(q(:))-1)>0.001)
        disp('q is not a probability distribution')
        return
    end
    
    
    for i=1:length(p(:))
        if(p(i)~=0)    
            D=D+p(i)*log(p(i)/q(i))/log(2);
        end
    end
end

