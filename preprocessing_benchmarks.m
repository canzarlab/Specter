addpath('LSC');
addpath('dimred');
% follow the original publication
% [data, rows, cols, entries] = mmread('/data/hoan/scRNAseq_clustering_benchmark/data/2Mcells_Trapnell/gene_count.txt');
% data = data';
% data: rows are samples, columns are genes

% data = csvread("/data/hoan/spectral_clustering_matlab/data/sim_big_data/SimDE001GuneqN10k_raw_count.csv");
% data = csvread("/data/hoan/spectral_clustering_matlab/data/sim_big_data/SimDE001GuneqN100k_raw_count.csv");
% data = csvread("/data/hoan/spectral_clustering_matlab/data/sim_big_data/SimDE001GuneqN200k_raw_count.csv");
% data = csvread("/data/hoan/spectral_clustering_matlab/data/sim_big_data/SimDE001GuneqN500k_raw_count.csv");
data = csvread("/data/hoan/spectral_clustering_matlab/data/sim_big_data/SimDE001GuneqN1m_raw_count.csv");

% n_select_genes = 2000; % number of selected genes in HV.
n_select_genes = 500; % number of selected genes in HV.
sample_size = 100000;
[m, n] = size(data)
subid = randperm(m);
fprintf('start here\n');
tic;
data = bsxfun(@rdivide,data, sum(data, 2));
if m <= sample_size
    [dump, varidx] = sort(var(data), 'descend'); % compute variance genes for subsample. 
else
    [dump, varidx] = sort(var(data(subid(1:sample_size),:)), 'descend'); % compute variance genes for subsample. 
end
topgenes = varidx(1:n_select_genes);
data = data(:,topgenes);
nA = 10000*bsxfun(@rdivide,data, sum(data,2));
nA = log(nA + 1);
% scale gene to have 0 mean and one variance
nA = zscore(nA);
% apply fast PCA
nA = fpca(nA, 100);
toc;

size(nA)
clear dump;
