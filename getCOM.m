function [COM]=getCOM(values)

	numer=0;
	for i=1:length(values)
		if(~isnan(values(i)))
			numer=numer+values(i)*i;
		end
	end

	denom=0;
	for i=1:length(values)
                if(~isnan(values(i)))
                        denom=denom+values(i);
                end
    end
        
    COM=numer/denom;
