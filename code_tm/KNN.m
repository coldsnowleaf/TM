function w = KNN(x, y, z, K)
% x: size k*m
% y: size k*n, database
% z: size l*n, database
% w: size l*m, output
% for each column in x, find its K nearest neighbors in the database y and
% output the average of the corresponding K columns in z


l_norm = 2;

m = size(x, 2);
n = size(y, 2);
k = size(x, 1);
l = size(z, 1);
if size(y, 1) ~= k || size(z, 2) ~= n
    error('Error in KNN(): the dimensions of the two inputs do not match.');
end
if K > n
    K = n;
end
    
w = zeros(l, m);
dist = zeros(m, n);
dist_Kmin_idx = zeros(K, m);
for i = 1:m
    for j = 1:n
        dist(i, j) = norm(x(:, i) - y(:, j), l_norm);
    end
    
    for t = 1:K
        dist_min = Inf;
        idx_min = 0;
        for j = 1:n
            if dist(i, j) < dist_min
                idx_min = j;
                dist_min = dist(i, j);
            end
        end
        dist_Kmin_idx(t, i) = idx_min;
        dist(i, idx_min) = Inf;
    end
    if K> 1
        w(:, i) = mean(z(:, dist_Kmin_idx(:, i))')';
    else
        w(:, i) = z(:, dist_Kmin_idx(:, i));
    end
end
