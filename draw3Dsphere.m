
numPts=100;
numPts=25;
numPts=20;
numPts=15;
x=1:numPts; y=1:numPts; z=1:numPts;

sphereRadSquared=3*(numPts/2)^2;
tol=100*(10/numPts);
tol=30*(10/numPts);
tol=50*(10/numPts);
c=zeros(length(x),length(y),length(z));
for xi=1:length(x)
    for yi=1:length(y)
        for zi=1:length(z)
            if(abs(xi^2+yi^2+zi^2-sphereRadSquared)<tol)
                c(xi,yi,zi)=1;
            end
        end
    end
end

plot4D(x,y,z,c,numPts)
