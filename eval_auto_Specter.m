function cls = eval_auto_Specter(fea, n_clusters, ensemble_size, mingamma)
    % Input: fea: expression data where rows are cells, collumns are PCs (genes)
    % n_clusters: number of clusters
    % ensemble_size: number of clusterings in the ensemble 
    % mingamma: minimum gaussion bandwidth (default: 0.1)
    
    [m, n] = size(fea);
    if m < 20000
        cls = eval_exact_Specter(fea, n_clusters, ensemble_size, mingamma);
    else
        n_neighbors = 5;
        cls = eval_fast_Specter(fea, n_clusters, ensemble_size, mingamma, n_neighbors);
    end

end
