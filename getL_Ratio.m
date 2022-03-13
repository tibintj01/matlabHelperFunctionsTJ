function [L_Ratio,L] = getL_Ratio(nonClusterD2,nCluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%assuming features are multivariate normally distributed, m-distance follows a chi-squared
%distribution with degrees freedom equal to the number of features 
%The CDF of this distribution gives the probability of finding an m-distance equal or lower
%than a given value
%The L ratio is the sum of the complement of this probability for each non-clustered point
%(probability of finding a clustered point greater than that distance), normalized by the
%number of spikes in the cluster 
%see Schmitzer-Torber et al. 2005 for further description
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-09-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFeatures=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing getL_Ratio fxn output.........')


L=0
tic
for i=1:length(nonClusterD2)
	L=L+(1-chi2cdf(nonClusterD2(i),numFeatures));
end

L_Ratio=L/nCluster;

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

