function cls = eval_fast_Specter(fea, n_clusters, ensemble_size, mingamma, n_neighbors)
    [m, n] = size(fea);
    samples_size = [n_clusters*100.0/sqrt(m)]; 
     % apply pre-processing
    params.HV = 1; % use HV seletion (0: no)
    params.PCA = 1;
    params.print = 1; % print results
    params.mode = 0; % 0: ensemble
    params.n_clusters = n_clusters;
    [opts, fea] = learn_LSC(fea, params);

    opts.kmMaxIter = 30;
    opts.maxIter = 100;
    if m > 500000
        opts.mode = 'random'
        opts.kmMaxIter = 10;
        opts.maxIter = 50;
    end
    N = ensemble_size; 
    fprintf('Ensemble size: %i\n', N);
    optsVec = repmat(opts, N);
    minsigma = mingamma; 

    clusters = zeros(N, m);
        
    parfor i=1:1:N
        optsVec(i).r = max(3, round(opts.r*(rand*0.5+0.8)));
        clusters(i,:) = LSC_eigen(fea, n_clusters, optsVec(i), rand*0.1 + minsigma); % rand*(1.5)+0.5 % VERY GOOD
    end

    n_samples = min(m, round(n_clusters*sqrt(m)));
 
    fprintf('run selective_sampling \n');
    cl = selective_sampling(clusters,n_clusters, n_samples); 

    selective_clusters = clusters(:, cl);
    unlabel_cl = setdiff(1:m, cl);
    trainData = fea(cl,:);
    baseCL = [cl unlabel_cl];
   
    %% --------------------- compute ensemble clustering -----------------------------------%% 
    ensemble = evalCOAL(selective_clusters, n_clusters); % run hierachical clustering of AL on co-association matrix
  
    fprintf('use kNN algorithm with k = %i \n', n_neighbors);
    model = fitcknn(trainData, ensemble,'NumNeighbors', n_neighbors); % nn = 5 
    testData = fea(unlabel_cl, :);
    prediction = predict(model,testData);
    ensemble = [ensemble; prediction];
    ensemble = sortrows([baseCL' ensemble], [1]);
    ensemble = ensemble(:,2);
    
    clear clusters baseCL cl optsVec prediction selective_clusters testData trainClass trainData unlabel_cl;
    cls = ensemble;

end
