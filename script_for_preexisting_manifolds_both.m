tot_num = 1;
P = 36;
avg_cns = zeros(P,tot_num);
avg_rads = zeros(P, tot_num);
avg_dims = zeros(P, tot_num);
avg_caps_theory = zeros(1,tot_num);
phi = zeros(P,P, tot_num);
lambda = zeros(P,P,tot_num);
mean_dots = zeros(P,tot_num);
exp_rest = zeros(1, tot_num);
thetas_rest = {};
tot_num_supp_vec = zeros(P, tot_num);


kappas = [0, 5, 25, 140, 529];

kappas_bw_mean = [16];
kappas_bw = [0,0.5, 2, 7.5, 15];
contrast = [0];
trials = 0:4;
count = 1;
for kappa = 1:length(kappas)
    for kappa_bw_mean = 1:length(kappas_bw_mean)
    for kappa_bw = 1:length(kappas_bw)
        for c=1:length(contrast)
            for t = trials
               [avg_cns(:,count), avg_rads(:,count), avg_dims(:,count), avg_caps_theory(count), phi( :,:,count), lambda(:,:,count), mean_dots(:,count),exp_rest(count), thetas_rest{count}, tot_num_supp_vec(:,count), corrs_per_dist(count,:)] = run_manifold_preexisting_both(kappas(kappa), kappas_bw(kappa_bw),kappas_bw_mean(kappa_bw_mean),contrast(c), t);
               count = count+1;
            end
        end
    end
    end
end

    