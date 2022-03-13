function []=compileAndRunMatlab(functionNamePath,inputCellMatPath,idPhrase)
%Tibin John
%Omar Ahmed Lab, Fall 2017

%intialize compiler function input, output
pathSplit=strsplit(functionNamePath,'/');
functionName=pathSplit{end};
functionName=functionName(1:(end-2));

functionDir=fullfile(pathSplit{1:(end-1)});

timestamp=getTimeStamp();
compiledDir=fullfile(pwd,['compiledDir_' timestamp '_' functionName '_' idPhrase]);

if(~isdir(compiledDir))
	mkdir(compiledDir)
end

recompile=1;

%compile matlab script
if(recompile)
	disp('compiling matlab script into bash executable........')
	tic
	mcc('-m','-R','-nodisplay','-R','-singleCompThread','-d',compiledDir,functionNamePath)
	toc
end

%#create gnu parallel cmds (each row of cell array is parallel instance, each column is parameter in argument list)
inputCellMat=load(inputCellMatPath);
inputCellArray=inputCellMat.inputCellArray;
cmdList=[];
count=0;
for parRowNum= 1:size(inputCellArray,1)
	
	cmd=sprintf('export MCR_CACHE_ROOT=/scratch/ojahmed_fluxm/tibintj/%d; %s %s', parRowNum, fullfile(compiledDir,['run_' functionName '.sh']), '/sw/sph/centos7/mcr/9.0.1')

	for argNum= 1:size(inputCellArray,2)
		typeChar=getTypeChar(inputCellArray{parRowNum,argNum});
		cmd=[cmd sprintf(sprintf(' %s',typeChar),inputCellArray{parRowNum,argNum})];
	end
	
	cmdList=[cmdList cmd sprintf('\n')];	
end
cmdList


%create pbs file using this executable from template including GNU parallel
cmdFileName=sprintf('%s_%s_gnuCmdsMatlab.txt',functionName,timestamp);
pbsFileName=sprintf('%s_%s_MatlabJobs.pbs',functionName,timestamp);

gnuCmdsFile=fopen(cmdFileName,'w');
fprintf(gnuCmdsFile,cmdList);
fclose(gnuCmdsFile)

pbsTemplateStr=fileread('/nfs/turbo/lsa-ojahmed/gnuParPbsTemplate.pbs')

%walltime='12:30:00'
walltime='02:30:00'
%walltime='08:40:00'

pbsStr=strrep(pbsTemplateStr,'JOBNAME',functionName);
pbsStr=strrep(pbsStr,'WALLTIME_STR',walltime);
pbsStr=strrep(pbsStr,'COMPILEDDIRPATH',fullfile(compiledDir));
pbsStr=strrep(pbsStr,'GNU_CMDFILEPATH',fullfile(pwd,cmdFileName))

newPBSFile=fopen(pbsFileName,'w');
fprintf(newPBSFile,pbsStr);
fclose(newPBSFile)

%submit job
jobID=deblank(evalc(sprintf('!qsub %s',pbsFileName)))

