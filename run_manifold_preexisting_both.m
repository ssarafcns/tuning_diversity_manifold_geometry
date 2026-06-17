function [ avg_cn, avg_r, avg_d, avg_ct, phi, lambda, mean_dots,exp_rest,thetas_rest, tot_num_supp_vec,corr_per_dist] =  run_manifold_preexisting_both(kappa, kappa_bw, kappa_bw_mean,contrast,trial)

kappa_string = num2str(kappa);

    kappabw_string = num2str(kappa_bw);
    if kappa_bw ==0.5
        kappabw_string='pt5';
    elseif kappa_bw == 7.5
        kappabw_string = '7pt5';
    end
contrast_string = num2str(100*contrast);
if contrast< 1
    contrast_string = ['pt', num2str(contrast*100)];
end

load(['het_',kappa_string,'_bw_', kappabw_string,'_mean_bw_',num2str(kappa_bw_mean),'_trial_',num2str(trial),'.mat']); % add contrast_pt0 to name string if contrast variations exp


[r, d, n, w, c_t, c_s, avg_cn, avg_r, avg_d, avg_ct, phi, lambda, mean_dots,exp_rest, thetas_rest,tot_num_supp_vec, corr_per_dist] = manifold_heterogeneity_FR_parallel_preexist_2(data_mat(:,:,1:36),1:36);
exp_rest = 0;
thetas_rest = 0;
tot_num_supp_vec = 0;
corr_per_dist = 0;
avg_ct = 0;
mean_dots = 0;
end