function [outputStruct] = plot2DGaussianContourGivenPts(pointCloud,dim1,dim2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-10-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lineWidth=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=pointCloud(:,[dim1 dim2]);
gm=fitgmdist(X,1);

mu = gm.mu;
Sigma = gm.Sigma;

%[.25 .3; .3 1];

%x1 = -3:.2:3; x2 = -3:.2:3;

currXlim=xlim;
currYlim=ylim;

%x1=linspace(currXlim(1),currXlim(2),1000);
%x2=linspace(currYlim(1),currYlim(2),1000);

%distXmin=prctile(X(:,1),1);
%distXmax=prctile(X(:,1),99);
%distYmin=prctile(X(:,2),1);
%distYmax=prctile(X(:,2),99);
distXmin=prctile(X(:,1),0.0000001);
distXmax=prctile(X(:,1),99.9999999);
distYmin=prctile(X(:,2),0.0000001);
distYmax=prctile(X(:,2),99.9999999);
%distXmin=prctile(X(:,1),0);
%distXmax=prctile(X(:,1),100);
%distYmin=prctile(X(:,2),0);
%distYmax=prctile(X(:,2),100);

x1=linspace(distXmin,distXmax,1000);
x2=linspace(distYmin,distYmax,1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing plot2DGaussianEllipse fxn output.........')

[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));
fMax=max(F(:))

F=F/sum(F(:)); %probability distribution or count?
normFMax=max(F(:))
%mvncdf([0 0],[1 1],mu,Sigma);
hold on
%contour(x1,x2,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
%contour(x1,x2,F,prctile(F(:),[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]*100));
%[C,h]=contour(x1,x2,F,prctile(F(:),[0.2:0.1:0.8]*100));
%[C,h]=contour(x1,x2,F,prctile(F(:),[0.25:0.25:0.75]*100));
%[C,h]=contour(x1,x2,F,prctile(F(:),[0.30:0.20:0.70]*100));
%[C,h]=contour(x1,x2,F,prctile(F(:),[0.5 0.6 0.8 ]*100));
[C,h]=contour(x1,x2,F,prctile(F(:),[0.7 0.7]*100),'b--','LineWidth',lineWidth);
hold on
[C,h]=contour(x1,x2,F,prctile(F(:),[0.8 0.8]*100),'b--','LineWidth',lineWidth);
[C,h]=contour(x1,x2,F,prctile(F(:),[0.9 0.9]*100),'b--','LineWidth',lineWidth);
%fds
%[C,h]=contour(x1,x2,F,prctile(F(:),[0.10:0.40:0.90]*100));
%[C,h]=contourf(x1,x2,F,prctile(F(:),[0.25:0.25:0.75]*100));

%[C,h]=contour(x1,x2,F,prctile(X(:),[0.4:0.1:0.7]*100));
%h.LineWidth=5;
%contour(x1,x2,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999],'LineWidth',5);
%contour(x1,x2,F,'LineWidth',5,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);

xlabel('x'); ylabel('y');
%saveas(gcf,'testNormalizedPDFcontour.tif')
%fds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

