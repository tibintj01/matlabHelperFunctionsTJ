function p=getGaussCurve(x,mu,sig)

    exponent=-0.5*((x-mu)/sig).^2;
    mag=1/(sig*sqrt(2*pi));
    
    p=mag*exp(exponent);