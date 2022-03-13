function [] = showFieldImage(currData,di) 
        if(di==2)
 
       
               if(iscell(currData.imagePaths))
                   imshow(currData.imagePaths{1})
               else
                   imshow(currData.imagePaths)
               end
 
        elseif(di==1)
  
                if(iscell(currData.imagePaths))
                   imshow(currData.imagePaths{2})
               else
                   imshow(currData.imagePaths)
                end

        end