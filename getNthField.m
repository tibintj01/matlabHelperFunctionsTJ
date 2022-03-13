function [nthFieldName,nthFieldVal]=getNthField(structVar,n)
                %structVar
		%n
		fieldNames=fieldnames(structVar);
                nthFieldName=fieldNames{n};
                nthFieldVal=structVar.(fieldNames{n});
