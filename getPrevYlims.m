function [prevYlims] =getPrevYlims()
	 prevYlims=ylim
	if(prevYlims(1)==0 && prevYlims(2)==1)
		
		prevYlims(1)=Inf;
		prevYlims(2)=-Inf;
	end
