function [global_mean_norm, dimension, norm_x0, gmsub_norm, raw_radius, radius, gmsub_radius, width, gmsub_width, capacity] = calculateManifoldDimensionRadius_Gaussian(points,global_mean, numP)
% points is MxN, where M is the number of trials, and N is the number
% of neurons

centroid = mean(points,1);
centered_points = points-repmat(centroid,length(points(:,1)),1);

global_mean_norm = norm(global_mean);
norm_x0 = norm(centroid);
gmsub_norm = norm(centroid - global_mean);


[U,S,V] = svd(centered_points);
eigs = sum(S,2).^2/(length(points(:,1))-1);
lambdas = sum(S,2);

%dimension = sum(lambdas)^2/sum(lambdas.^2);  % participation ratio
%raw_radius = sqrt(sum(lambdas.^2));%mean(lambdas);
%dimension = (sum(lambdas.^2)^2) / sum(lambdas.^4); %sorscher PR without 1/N
%
[raw_radius, dimension] = compute_gauss_RD(centered_points', 1000);

radius = raw_radius/norm_x0;
gmsub_radius = raw_radius/gmsub_norm;


width = radius*sqrt(dimension);
gmsub_width = gmsub_radius*sqrt(dimension);

capacity = (1+ gmsub_radius^(-2))/dimension;

end