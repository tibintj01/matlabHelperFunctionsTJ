function [finalCmap] = getBiColorMap(n)
%get bicolormap

%data=load('Figure4b_cmap1.mat');
data=load('Figure4b_cmap2.mat');

numColors=size(data.mycmap,1);

  mycmap=data.mycmap;
  
  
  

  
if(~exist('n','var'))
  finalCmap=mycmap;

else
    finalCmap=NaN(n,3);
      newColorRows=round(linspace(1,numColors,n));
      
    for i=1:length(newColorRows)
        finalCmap(i,:)=mycmap(newColorRows(i),:);
        %colorCount=colorCount+1;
    end
    
    disp('')
    
end



