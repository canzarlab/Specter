clear all;
addpath('dimred');
addpath('LSC');
addpath('utils');
format short g;

% read mRNA data
fea = csvread('data/cbmc_rna_pca.csv'); % read a csv file
% read ADT data
fea_adt = csvread('data/cbmc_adt.csv'); % read a csv file

[m, n] = size(fea);
n_clusters = 15;
% get the parameter for model 
params.HV = 1; % use HV seletion (0: no)
params.PCA = 1;
params.print = 1; % print results
params.mode = 0; % 0: ensemble
params.n_clusters = n_clusters;
[opts, fea] = learn_LSC(fea, params);

% compute ensemble
N = 200; 
M = 200;
fprintf('Ensemble size: %i\n', N + M);
clusters = zeros(N + M, m);

% for RNA count
opts_rna = opts;
optsVec_rna = repmat(opts_rna, N);
parfor i=1:1:N
    optsVec(i).r = max(3, round(opts.r*(rand*0.5+0.8)));
    clusters(i,:) = LSC_eigen(fea, n_clusters, optsVec_rna(i), rand*1.1 + 1.4); 
end
% for ADT count 
optsVec = repmat(opts, M);
parfor i = N+1:1:(N+M)
    optsVec(i-N).r = max(3, round(opts.r*(rand*0.5+0.8)));
    clusters(i,:) = LSC_eigen(fea_adt, n_clusters, optsVec(i-N), rand*0.1 + 0.3); 
end

ensemble = evalCOAL(clusters, n_clusters); % run hierachical clustering of AL on co-association matrix

% save label 
writematrix(ensemble, strcat("output/cbmc_SCE_rna_labels.csv"));



