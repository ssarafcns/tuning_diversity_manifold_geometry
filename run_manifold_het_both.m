kappas = [0,5,25,140,529];
kappas_bw = [0, 0.5, 2, 7.5, 15];
contrast = [0];
kappas_bw_mean = [16];
Ns = [300];
num_trials = 50;
total_num_runs = length(Ns)*length(kappas)*length(kappas_bw)*num_trials;
capacity_theory = zeros(1,total_num_runs);
capacity_sim = zeros(1,total_num_runs);
raw_rads = zeros(1,total_num_runs);
dims = zeros(1,total_num_runs);
norms = zeros(1,total_num_runs);
widths = zeros(1,total_num_runs);
avg_cns = zeros(1, total_num_runs);
avg_rads = zeros(1,total_num_runs);
avg_dims = zeros(1, total_num_runs);
avg_caps_theory = zeros(1, total_num_runs);
means_amp = zeros(1,total_num_runs);
means_bw = zeros(1, total_num_runs);
vars_amp = zeros(1, total_num_runs);
vars_bw = zeros(1, total_num_runs);
lambda = zeros(36,36,total_num_runs);


count = 1;
for kappa = 1:length(kappas)
for kappa_bw_mean = 1:length(kappas_bw_mean)
    for kappa_bw = 1:length(kappas_bw)
        for c=1:length(contrast)
        for trial=0:num_trials-1
    plot_tc = 0;
    if trial ==0
        plot_tc = 1;
    end

    %[r, d, n, w, c_t, c_s,a_cn, a_r, a_d, a_ct,var_amp, mean_amp, var_bw,mean_bw] = manifold_heterogeneity_FR_parallel_both_calc_03_11_24(kappas(kappa),kappas_bw(kappa_bw),kappas_bw_mean(kappa_bw_mean),contrast(c), Ns(1), trial,0); %% UNCOMMENT FOR CONTRAST VARIATIONS
    [r, d, n, w, c_t, c_s,a_cn, a_r, a_d, a_ct, var_amp, mean_amp, var_bw, mean_bw,lambda(:,:,count)] = manifold_heterogeneity_FR_parallel_both_calc_01_05_24(kappas(kappa),kappas_bw(kappa_bw),kappas_bw_mean(kappa_bw_mean), Ns(1), trial,0);
    capacity_theory(count) = a_ct;
    avg_caps_theory(count) = a_ct;
    capacity_sim(count) = c_s;
    raw_rads(count) = r;
    dims(count) =  d;
    norms(count) = n;
    widths(count) = w;
    avg_cns(count) = a_cn;
    avg_rads(count) = a_r;
    avg_dims(count) = a_d;
    means_amp(count) = mean_amp;
    vars_amp(count) = var_amp;
    means_bw(count) = mean_bw;
    vars_bw(count) = var_bw;
    count = count+1;
        end
        end
    end
end
end
