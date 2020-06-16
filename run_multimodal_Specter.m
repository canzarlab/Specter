clear all;
addpath('dimred');
addpath('LSC');
addpath('utils');
format short g;

% read mRNA data
fea = csvread('data/pbmc_process_RNA_pca.csv'); % read a csv file
% read ADT data
fea_adt = csvread('data/pbmc_process_ADT.csv'); % read a csv file

[m, n] = size(fea);
gamma_rna = 0.4; 
gamma_adt = 0.3;
n_clusters_rna = 16; 
n_clusters_adt = 16;
% get the parameter for model 
params.HV = 1; % use HV seletion (0: no)
params.PCA = 1;
params.print =1; % print results
params.mode = 2; 
params.n_clusters = n_clusters_rna;
[opts, fea] = learn_LSC(fea, params);

% compute ensemble
N = 200; 
M = 200;
fprintf('Ensemble size: %i\n', N + M);
clusters = zeros(N + M, m);

% for RNA count
tic;
opts_rna = opts;
optsVec_rna = repmat(opts_rna, N);
parfor i=1:1:N
    optsVec(i).r = max(3, round(opts.r*(rand*0.5+0.8)));
    clusters(i,:) = LSC_eigen(fea, n_clusters_rna, optsVec_rna(i), rand*1.1 + gamma_rna); 
end
% for ADT count 
optsVec = repmat(opts, M);
parfor i = N+1:1:(N+M)
    optsVec(i-N).r = max(3, round(opts.r*(rand*0.5+0.8)));
    clusters(i,:) = LSC_eigen(fea_adt, n_clusters_adt, optsVec(i-N), rand*0.1 + gamma_adt); 
end

ensemble = evalCOAL(clusters, 16); % run hierachical clustering of AL on co-association matrix
toc;

% save label 
writematrix(ensemble, strcat("output/pbmc_results.csv"));



