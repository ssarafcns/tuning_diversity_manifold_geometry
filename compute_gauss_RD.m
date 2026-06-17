function [Rg, Dg] = compute_gauss_RD(S, n_t)
  
%     Additional function to calculate Gaussian radius and dimension of a manifold. 
%     NOTE: S should be centered on the origin.
%     Args:
%         S: Manifold matrix 
%         n_t: Number of gaussian vectors to sample per manifold
%         
%     Returns:
%         Rg: Gaussian Radius of the input manifold
%         Dg: Gaussian Dimension of the input manifold
  
    stream = RandStream('mt19937ar', 'Seed', 1000);
    [N, m] = size(S);
    t_vec = randn(stream,N, n_t);
    smax = zeros(N,n_t);
    smax_norm_all = [];
    tShat_all = [];
    for i =1:n_t
        t = t_vec(:,i);  %nx1

        %print("t",t.shape)
        A = t'*S; %(1xn)(nXm)
        [tSmax, imax] = max(A);

        smax = S(:,imax);
        smax_norm_all = [smax_norm_all, norm(smax,2)];
        tShat_all = [tShat_all,  t'*(smax/norm(smax,2)) ];
    end
    Rg_sq = mean(smax_norm_all.^2);
    Rg = sqrt(Rg_sq);
    Dg = mean(tShat_all.^2);
   end