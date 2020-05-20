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
adt_K = 15;
rna_K = 16;
n_clusters = 18;
% get the parameter for model 
params.HV = 1; % use HV seletion (0: no)
params.PCA = 1;
params.print = 1; % print results
params.mode = 0; % 0: ensemble
params.n_clusters = n_clusters;
[opts, fea] = learn_LSC(fea, params);

% compute ensemble
N = 50; 
M = 20;
fprintf('Ensemble size: %i\n', N + M);
clusters = zeros(N + M, m);

% for RNA count
tic;
opts_rna = opts;
optsVec_rna = repmat(opts_rna, N);
parfor i=1:1:N
    % optsVec(i).r = max(3, round(opts.r*(rand*1.5+10.0)));
    % optsVec(i).r = 10;
    % clusters(i,:) = LSC_eigen(fea, rna_K, optsVec_rna(i), rand*1.1 + 1.4); 
    % clusters(i,:) = LSC_eigen(fea, rna_K, optsVec_rna(i), rand*0.5 + 0.9); 
    % clusters(i,:) = LSC_eigen(fea, rna_K, optsVec_rna(i), rand*0.5 + 0.8); %ok for k = 18
    clusters(i,:) = LSC_eigen(fea, rna_K, optsVec_rna(i), rand*0.5+ 0.9); 
end
% for ADT count 
optsVec = repmat(opts, M);
parfor i = N+1:1:(N+M)
    % optsVec(i-N).r = max(3, round(opts.r*(rand*0.5+0.8)));
    % clusters(i,:) = LSC_eigen(fea_adt, adt_K, optsVec(i-N), rand*0.5 + 0.3);  ok
    % clusters(i,:) = LSC_eigen(fea_adt, adt_K, optsVec(i-N), rand*0.1 + 0.1); % find Eryth but not DCs
    clusters(i,:) = LSC_eigen(fea_adt, adt_K, optsVec(i-N), rand*0.1 + 0.1); 
end

ensemble = evalCOAL(clusters, n_clusters); % run hierachical clustering of AL on co-association matrix
toc;

% save label 
writematrix(ensemble, strcat("output/cbmc_specter_adtK_", num2str(adt_K), "_rnaK_", num2str(rna_K), "_labels_v4.csv"));
% writematrix(ensemble, strcat("output/cbmc_specter_rna.csv"));



