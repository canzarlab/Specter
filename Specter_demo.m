clear all;
% load all necessary library
addpath('dimred');
addpath('LSC');
addpath('utils');
format short g;

%% read gene expression data: rows are cells, collumns are features (PCs, genes)
data = csvread("data/zeisel_pca.csv");  % e.g., read zeisel data (reduced data by PCA)

%% choose parameters for algorithm
n_clusters = 7; % number of clusters 
ensemble_size = 200; % ensemble sizes (default: 200)
mingamma = 0.1; % minimum gaussian bandwidth (default: 0.1)

%% We have developed two version of Specter. Let Specter to decide the algorithm
tic;
specter_labels = eval_auto_Specter(data, n_clusters, ensemble_size, mingamma);
toc;

%% Evaluate the solution (ARI score)
labels_true = csvread("data/zeisel_pca_labels.csv");
score = eval_rand(labels_true, specter_labels);


%% To run a specific an algorithm of Specter:
% exact algorithm with time complexity of O(n^2) where n is number of cells, or
% fast Specter with time complexity of O(n)

% exact Specter
exact_specter_labels = eval_exact_Specter(data, n_clusters, ensemble_size, mingamma);
score_v2 = eval_rand(labels_true, exact_specter_labels);

% fast Specter
n_neighbors = 5; %% parameter for k-nearest neighbor algorithm (default: 5)
fast_specter_labels = eval_fast_Specter(data, n_clusters, ensemble_size, mingamma, n_neighbors);
score_v3 = eval_rand(labels_true, exact_specter_labels);
