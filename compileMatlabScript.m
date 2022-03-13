function []=compileAndRunMatlab(functionNamePath)

pathSplit=strsplit(functionNamePath,'/');
functionName=pathSplit{end};
functionName=functionName(1:(end-2));

functionDir=fullfile(pathSplit{1:(end-1)});

compiledDir=fullfile(['compiledDir_' functionName]);

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

cmdList=[];

count=0;
for feature1 = 1:18

	for feature2= 1:18
		count=count+1;
		cmd=sprintf('export MCR_CACHE_ROOT=/scratch/ojahmed_fluxm/tibintj/%d; %s %s', count, fullfile(pwd,compiledDir,['run_' functionName '.sh']), '/sw/sph/centos7/mcr/9.0.1')

		cmd=[cmd sprintf(' %d %d',feature1,feature2)];

		cmdList=[cmdList cmd sprintf('\n')];	
	end
end
cmdList

gnuCmdsFile=fopen('gnuCmdsMatlab.txt','w');
fprintf(gnuCmdsFile,cmdList);
fclose(gnuCmdsFile);
