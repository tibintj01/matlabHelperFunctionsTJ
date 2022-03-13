function [concatFileHandle]=loadConcatFile(chanNum,concatDir)

	if(~exist('concatDir'))
                %disp('Using directory for MG49 session 3......')
                concatDir='C:\Users\tibin\spikeDynamicsAnalysisTibin\session3Mat';
        end

        filePath=fullfile(concatDir,sprintf('concatChan%d',chanNum), 'lfpOriginalConcat.mat');


        concatFileHandle=matfile(filePath);

