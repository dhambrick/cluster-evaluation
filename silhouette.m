1;
function [b] = CalcIntraClusterDistance(clusters,cluster_idx,data_idx,distances)
               n_clusters = numel(unique(clusters));
               cluster = find(clusters == cluster_idx);
               cluster_size = length(cluster);
               dist_i = distances(data_idx,cluster);
               b = mean(dist_i);
end

function [silhouettes] = CalcSilhouettes(data,clusters)
        n_of_data = size(data);
        uniq_clusters = unique(clusters);
        n_clusters = numel(uniq_clusters);
        a = zeros(n_of_data(1),1);
        b = zeros(n_of_data(1),1);
        silhouettes = zeros(n_of_data(1),1);
        if(n_clusters== 1)
            error('Data needs k > 1 clusters in order to calculate silhouettes')
        end
          
        %Sort the cluster vector 
        [sorted_clusters, sorted_idx] = sort(clusters);
        %Calculate all the squared euclidean distances
        distances = zeros(n_of_data(1),n_of_data(1));
          
        % We want to calculate the dist from each point to every other point in data
        % Normally this would be O(N^2) , but since we exclude dist(d_i,d_i) and since the
        % distances are symmetric, this really reduces down to O(N^2/2)
          
        [r_idx , c_idx] = ind2sub([n_of_data,n_of_data],find(tril(not(eye(n_of_data(1))))==1));
         
        for idx = 1:length(r_idx)
            r = r_idx(idx);
            c = c_idx(idx);
            distances(r , c) = norm(data(r,:)-data(c,:),2);
            distances(c , r) = distances(r,c);
        endfor
          
        %For each data point, calculate average inter cluster distance:
        for idx = 1:n_of_data(1)
            rs_cluster = clusters(idx);
            % find all of the data points that share the same cluster as r
            rs_neighs = find(clusters == rs_cluster);
            n_neighs = numel(rs_neighs);
            %for the average simillarity
            distances(idx,rs_neighs);
            a(idx) = (1/(n_neighs-1))*sum(distances(idx,rs_neighs));
            other_clusters = uniq_clusters(uniq_clusters != rs_cluster);
            temp_b = zeros(length(other_clusters),1);
            for jdx = 1:length(other_clusters)
                temp_b(jdx) = CalcIntraClusterDistance(clusters,other_clusters(jdx),idx,distances);
            endfor
            b(idx) = min(temp_b);
        endfor
        for kdx = 1:length(a)
            n = b(kdx) - a(kdx);
            d = max(a(kdx),b(kdx));
            s(kdx) = n/d;
        endfor
        silhouettes = s;
endfunction

function [s_by_cluster , partioned_cluster_idx] = PartitionSilhouette(silhouettes,clusters)
    uniq_clusters = unique(clusters);
    n_uniq_clusters = numel(uniq_clusters);
    partioned_cluster_idx = {};
    s_by_cluster = {};
    for idx = 1:n_uniq_clusters
        cluster_id = uniq_clusters(idx);
        temp = {cluster_id,find(clusters == cluster_id)};
        partioned_cluster_idx(idx,:) = temp;
        s_by_cluster(idx) = silhouettes(cell2mat(temp(1,2)));
    endfor
endfunction

function GraphSilhouetteClusters(data_labels,silhouettes,clusters)
    [sil_by_clust,part_clust_idx] = PartitionSilhouette(silhouettes,clusters);
    n_clusters = length(part_clust_idx);
    clf
    for idx = 1:n_clusters
        subplot(n_clusters,1,idx)
        barh(cell2mat(sil_by_clust(idx)),"grouped")
    endfor
    
endfunction
