function []=compileAndRunMatlab(functionNamePath,inputCellMatPath,idPhrase)

pathSplit=strsplit(functionNamePath,'/');
functionName=pathSplit{end};
functionName=functionName(1:(end-2));

functionDir=fullfile(pathSplit{1:(end-1)});

timestamp=getTimeStamp();
compiledDir=fullfile(pwd,['compiledDir_' timestamp '_' functionName '_' idPhrase]);
%compiledDir=fullfile(['compiledDir_' functionName]);

if(~isdir(compiledDir))
	mkdir(compiledDir)
end

recompile=1;
if(recompile)
	disp('compiling matlab script into bash executable........')
	tic
	mcc('-m','-R','-nodisplay','-R','-singleCompThread','-d',compiledDir,functionNamePath)
	toc
end

%#create gnu parallel cmds

%#export MCR_CACHE_ROOT=/lscratch/$SLURM_JOBID; shell_scriptsgetInspActivity/run_getInspActivity.sh /usr/local/matlab-compiler/v91 171019-2_00Ctrl
%#export MCR_CACHE_ROOT=/scratch/ojahmed_fluxm/tibintj/; shell_scriptsgetInspActivity/run_getInspActivity.sh /usr/local/matlab-compiler/v91 171019-2_00Ctrl
%#/scratch/ojahmed_fluxm/tibintj/
%#this works: ./run_plotAllFeatureCombs.sh /sw/sph/centos7/mcr/9.0.1

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


cmdFileName=sprintf('%s_%s_gnuCmdsMatlab.txt',functionName,timestamp);
pbsFileName=sprintf('%s_%s_MatlabJobs.pbs',functionName,timestamp);

gnuCmdsFile=fopen(cmdFileName,'w');
fprintf(gnuCmdsFile,cmdList);
fclose(gnuCmdsFile)

pbsTemplateStr=fileread('/nfs/turbo/lsa-ojahmed/gnuParPbsTemplate.pbs')

walltime='12:00:00'

pbsStr=strrep(pbsTemplateStr,'JOBNAME',functionName);
pbsStr=strrep(pbsStr,'WALLTIME_STR',walltime);
pbsStr=strrep(pbsStr,'COMPILEDDIRPATH',fullfile(compiledDir));
pbsStr=strrep(pbsStr,'GNU_CMDFILEPATH',fullfile(pwd,cmdFileName))

newPBSFile=fopen(pbsFileName,'w');
fprintf(newPBSFile,pbsStr);
fclose(newPBSFile)

jobID=deblank(evalc(sprintf('!qsub %s',pbsFileName)))

