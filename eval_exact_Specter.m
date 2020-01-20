function cls = eval_exact_Specter(fea, n_clusters, ensemble_size, mingamma)
    [m, n] = size(fea);
    % apply pre-processing
    params.mode =0;
    params.HV = 1; % use HV seletion (0: no)
    params.PCA = 1;
    params.print = 1; % print results
    params.n_clusters = n_clusters;
    [opts, fea] = learn_LSC(fea, params);
    
    N = ensemble_size;
    
    clusters = zeros(N,m);
    optsVec = repmat(opts, N);
    parfor i=1:1:N
        optsVec(i).r = max(3, round(opts.r*(rand*0.5+0.8)));
        clusters(i,:) = LSC_eigen(fea, n_clusters, optsVec(i), rand*0.2 + mingamma); %  very good for co-association clustering
    end

    cls = evalCOAL(clusters, n_clusters); % run hierachical clustering of AL on co-association matrix

end
