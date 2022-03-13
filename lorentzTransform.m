function [xprime,tprime]=lorentzTransform(x,t,v,c)

    gamma=1/sqrt(1-(v^2/c^2));
    
    xprime=gamma*(x-v*t);
    tprime=gamma*(t-x*v/c^2);
end

