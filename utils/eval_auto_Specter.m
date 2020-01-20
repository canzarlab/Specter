function cls = eval_auto_Specter(fea, n_clusters, ensemble_size, mingamma)
    [m, n] = size(fea);
    if m < 20000
        cls = eval_exact_Specter(fea, n_clusters, ensemble_size, mingamma);
    else
        n_neighbors = 5;
        cls = eval_fast_Specter(fea, n_clusters, ensemble_size, mingamma, n_neighbors);
    end

end
