
function centroids = initCentroids(X, K)
    centroids = zeros(K,size(X,2)); 
    randidx = randperm(size(X,1));
    centroids = X(randidx(1:K), :);
  end

function indices = getClosestCentroids(X, centroids)
  K = size(centroids, 1);
  indices = zeros(size(X,1), 1);
  m = size(X,1);

  for i=1:m
    k = 1;
    min_dist = sum((X(i,:) - centroids(1,:)) .^ 2);
    for j=2:K
        dist = sum((X(i,:) - centroids(j,:)) .^ 2);
        if(dist < min_dist)
          min_dist = dist;
          k = j;
        end
    end
    indices(i) = k;
  end
end

function centroids = computeCentroids(X, idx, K)

  [m n] = size(X);
  centroids = zeros(K, n);
  
  for i=1:K
    xi = X(idx==i,:);
    ck = size(xi,1);
    # centroids(i, :) = (1/ck) * sum(xi);
    centroids(i, :) = (1/ck) * [sum(xi(:,1)) sum(xi(:,2))];
  end
end

function  [assignments, centers] = kmeans(X, k, centers = 0, maxiter = 200)
	if (centers == 0)
		centerRows = randperm(size(X)(1));
		centers = X(centerRows(1:k), :);
	endif
	numOfRows = length(X(:,1));
	numOfFeatures = length(X(1,:));
	assignments = ones(1, numOfRows);

	for iter = 1:maxiter
		clusterTotals = zeros(k, numOfFeatures);
		clusterSizes = zeros(k, 1);
		for rowIx = 1:numOfRows
			minDist = realmax;
			assignTo = 0;
			for centerIx = 1:k 
				% Euclidian distance is used.
				dist = sqrt(sum((X(rowIx, : ) - centers(centerIx, :)).^2));
				if dist < minDist
					minDist = dist;
					assignTo = centerIx;
				endif
			endfor
			assignments(rowIx) = assignTo;

			% Keep these information to calculate cluster centers.
			clusterTotals(assignTo, :) += X(rowIx, :);
			clusterSizes(assignTo)++;
		endfor

		% This process is called 'singleton' in terms of Matlab. 
		% If a cluster is empty choose a random data point as new 
		% cluster cener.
		for clusterIx = 1:k
			if (clusterSizes(clusterIx) == 0)
				randomRow = round(1 + rand() * (numOfRows - 1) );
				clusterTotals(clusterIx, :) =  X(randomRow, :);
				clusterSizes(clusterIx) = 1;
			endif
		endfor

		newCenters = zeros(k, numOfFeatures);
		for centerIx = 1:k 
			newCenters(centerIx, :) = clusterTotals(centerIx, : ) / clusterSizes(centerIx);
		endfor
	
		diff = sum(sum(abs(newCenters - centers)));
	
		if diff < eps
			%disp('Centers are same, which means we converged before maxiteration count. This is a good thing!')
			break;
		endif
	
		centers = newCenters;
	endfor	
	assignments = assignments';
	%printf('iter: %d, diff: %f\n', iter, diff);
endfunction

X = magic(5)

[a,c] = kmeans(X,2,centers=0,maxiters=200)

find(a==1)

function [b] = CalcIntraClusterDistance(clusters,cluster_idx,data_idx,distances)
               n_clusters = numel(unique(clusters))
               cluster = find(clusters == cluster_idx)
               cluster_size = length(cluster)
               dist_i = distances(data_idx,cluster)
               b = mean(dist_i)
               

function [silhouettes] = CalcSilhouettes(data,clusters)
          n_of_data = size(clusters)
          n_clusters = numel(unique(clusters))
          a = zeros(n_of_data)
          
          silhouettes = zeros(n_of_data);
          if(n_clusters)== 1);
              error('Data needs k > 1 clusters in order to calculate silhouettes')
          endif
          
          %Sort the cluster vector 
          [sorted_clusters, sorted_idx] = sort(clusters);
          %Calculate all the squared euclidean distances
          distances = zeros(n_of_data(1),n_of_data(1));
          
          % We want to calculate the dist from each point to every other point in data
          % Normally this would be O(N^2) , but since we exclude dist(d_i,d_i) and since the
          % distances are symmetric, this really reduces down to O(N^2/2)
          
          [r_idx , c_idx] = ind2sub([n_of_data,n_of_data],find(tril(not(eye(n_of_data(1))))==1));
          disp(r_idx)
          for idx = 1:length(r_idx)
              r = r_idx(idx);
              c = c_idx(idx);
              distances(r , c) = norm(data(r),data(c),2).^2;
              distances(c , r) = distances(r,c);
          endfor
          disp(distances)
          %For each data point, calculate average inter cluster distance:
          for idx = 1:n_of_data(1)
              rs_cluster = clusters(idx);
              % find all of the data points that share the same cluster as r
              rs_neighs = find(clusters == rs_cluster)
              n_neighs = numel(rs_neighs)
              %for the average simillarity
              a(idx) = (1/(n_neighs-1))*sum(distances(idx,rs_neighs))
              
          endfor
          
          
endfunction
          
          

s = CalcSilhouettes(X,[1])


