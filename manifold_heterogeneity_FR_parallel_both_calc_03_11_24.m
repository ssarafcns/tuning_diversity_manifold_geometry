function [raw_rad, dimension, gmsub_norm, width, capacity_theory, capacity_sim, avg_cns,avg_rs, avg_ds, avg_caps_theory, var_amp, mean_amp, var_bw, mean_bw] = manifold_heterogeneity_FR_parallel_both_calc_03_11_24(kappa, kappa_bw, kappa_bw_mean,contrast, N, trial, plot_tcs)

rho = 1;
c = 0.4;
a = 40;
sig = 1/sqrt(2);
fmax = 60;
fref = 3.3;
rng(trial+105) 
factor = 10/10;
leftover = N*factor;
k = 1:1:leftover;
%make preferred angles
phi = -1*pi*(leftover+1)/leftover + 2*pi/leftover*k;
phi = mod(phi, 2*pi);
phi_2 = deg2rad(180)*ones(1,N-leftover);
rand_phi_2 = deg2rad(rand(1,N-leftover)*20 - 10);
phi_2 = phi_2 + rand_phi_2;
phi = [phi, phi_2];
phi = mod(phi, 2*pi);

% make covariance matrix

Cor = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i ==j
           Cor(i,j) = 1;
        else
            diff_theta_normal = abs(phi(i)-phi(j));
            diff_theta_across = 2*pi-(abs(phi(i)-phi(j)));
            diff_theta = min(diff_theta_normal, diff_theta_across);
            Cor(i,j) = c*exp( -1/rho * abs(diff_theta) );
        end
    end
end

start_theta = 0;
step_deg = 5;
end_theta = 175;
num_trials=50;
thetas = start_theta:step_deg:end_theta;


mu_amp = 30;
if kappa <-1  
    rand_var_p2t = ones(1,N)*kappa;
 else
    mean_amp = 0;
    rand_var_p2t = random('Lognormal',log(mu_amp^2/sqrt(kappa+mu_amp^2)),sqrt(log(kappa/(mu_amp^2)+1)),[1,N]);
    mean_amp = mean(rand_var_p2t);
end
var_amp = std(rand_var_p2t);
mean_amp = mean(rand_var_p2t);


x = acos(1+ log(1/sqrt(2))/2 );
if kappa_bw==0
  rand_var_gamma = x/deg2rad(kappa_bw_mean);
  mean_bw = kappa_bw_mean;
  var_bw = 0;
else
  new_scale = kappa_bw;
  k = kappa_bw_mean/new_scale;
  mean_bw = 0;
  rand_var_gamma = gamrnd(k,new_scale,1,N);
  mean_bw = mean(rand_var_gamma);
  mean_bw= mean(rand_var_gamma);
  var_bw = std(rand_var_gamma); 
  rand_var_gamma = x./deg2rad(rand_var_gamma);

end

bw_fac = rand_var_gamma;

%DS levels : 2=DS, 1= DS Bias, 0 = OS
ds_levels = rand(1,N);
ds_levels(ds_levels <.1) = 0;
ds_levels(ds_levels < .3) = 0;
ds_levels(ds_levels<1) = 0;

opp_phi = mod(phi +pi,2*pi);

hold_thetas = repmat(deg2rad(0:1:359)',1,N);
hold_thetas_w_fac = hold_thetas.*bw_fac;
f_0_hold = exp( (cos(hold_thetas_w_fac) - 1) / sig^2 );
f_0_hold_opp = exp( (cos(abs(pi-hold_thetas).*bw_fac) - 1) / sig^2 );
mat_max = f_0_hold < .02;
[~,inds_max] = max(mat_max, [], 1);
thetas_max = hold_thetas(inds_max);

counts = zeros(num_trials*length(thetas), N);
centroids = zeros(length(thetas),N);
ms = zeros(N, length(thetas));
for theta=thetas
  theta_rad = deg2rad(theta);
  diff_theta_normal = abs(theta_rad-phi);
  diff_theta_other = 2*pi - abs(theta_rad-phi);

  diff_theta = min(diff_theta_normal,diff_theta_other);

  diff_theta_opp_normal = abs(theta_rad-opp_phi);
  diff_theta_opp_other = 2*pi - abs(theta_rad-opp_phi);
  diff_theta_opp = min(diff_theta_opp_normal,diff_theta_opp_other);

  inds_over_max = diff_theta >=thetas_max;
  inds_over_max_opp = diff_theta_opp >= thetas_max;
  diff_theta = diff_theta.*bw_fac;
  diff_theta_opp = diff_theta_opp.*bw_fac;

  f_0_pref = rand_var_p2t.*exp( (cos(diff_theta) - 1) / sig^2 );
  f_0_pref(inds_over_max) =rand_var_p2t(inds_over_max)*.02;
  f_0_opp = zeros(1,N);
  f_0_opp(ds_levels==0) = rand_var_p2t(ds_levels==0).*exp( (cos(diff_theta_opp(ds_levels==0)) - 1) / sig^2 );
  f_0_opp(ds_levels==1) = 0.6*rand_var_p2t(ds_levels==1).*exp( (cos(diff_theta_opp(ds_levels==1)) - 1) / sig^2 );
  f_0_opp(inds_over_max_opp) = rand_var_p2t(inds_over_max_opp)*0.02;
  
  f_0 = f_0_pref;
  f_0(ds_levels<2) = max(f_0_pref(ds_levels<2), f_0_opp(ds_levels<2));

  hold_cts = zeros(num_trials, N);
  for n_t = 1:num_trials
    ctrst = 0.5*rand()+.55;
    m_0 = (ctrst^3/(ctrst^3+0.5))*f_0 + fref;
    ms(:,theta/step_deg+1) = m_0';

    sqrt_m_0 = real(sqrt(m_0));
    C = bsxfun(@times, sqrt_m_0,Cor);
    C = bsxfun(@times, sqrt_m_0', C);
  
    hold_cts((n_t-1)*1+1:n_t*1,:) = mvnrnd(m_0,C, 1);
  end
   theta_num =floor((theta-start_theta)/step_deg);
   counts(theta_num*num_trials+1:theta_num*num_trials+num_trials,:) = hold_cts;
   centroids(theta_num+1,:) = mean(counts(theta_num*num_trials+1:theta_num*num_trials+num_trials,:),1);
   t2tvar(theta_num+1,:) = var(counts(theta_num*num_trials+1:theta_num*num_trials+num_trials,:),1);
end
kappa_string = num2str(kappa);

kappabw_string = num2str(kappa_bw);
if kappa_bw == 0.5
    kappabw_string = ['pt', num2str(kappa_bw*10)];
end

if contrast <1
    contrast_string = ['pt', num2str(contrast*100)];
else
    contrast_string = num2str(contrast*100);
end
data_mat = reshape(counts', N, num_trials, length(thetas));
save(['het_',kappa_string,'_bw_', kappabw_string, '_mean_bw_', num2str(kappa_bw_mean),'_contrast_', contrast_string,'_trial_',num2str(trial),'.mat'], "data_mat")
num_thetas = length(thetas);
global_mean = mean(counts,1);

if plot_tcs
    fig = figure;
    for c=1:N
      plot(thetas, ms(c,:),'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5); hold on;
    end
    box off;
    ax = gca; 
    ylabel('response (spikes/s)')
    xlabel('direction (degrees)')
    xlim([-5,365])
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;    
    title(['kappa = ', num2str(kappa), ' kappa bw = ', num2str(kappa_bw)]);
    saveas(fig,['kappa = pt', num2str(kappa*100), ' kappa bw = ', num2str(kappa_bw)], 'pdf' )
   
    
    figure;
    bws = rad2deg(acos(0.5*log(1/sqrt(2)) +1) ./ bw_fac);
    hist(bws);
    max_c = max(centroids,[],1);
    min_c = min(centroids,[],1);
    title(['kappa = ', num2str(kappa), ' kappa bw = ', num2str(kappa_bw), ' mean amp', num2str(mean(max_c-min_c))]);
    

    figure;
    hist(rand_var_p2t);
    title(['kappa = ', num2str(kappa), ' kappa bw = ', num2str(kappa_bw)]);
    
    
end

avg_cns = 0;
avg_rs = 0;
avg_ds = 0;
avg_caps_theory = 0;

num_to_get_geom = length(thetas);
for t=1:num_to_get_geom
[~,dimension, ~,gmsub_norm,raw_rad,~,~, ~, width, capacity_theory] = calculateManifoldDimensionRadius_Gaussian(counts((t-1)*num_trials+1:t*num_trials,:),mean(counts,1), num_to_get_geom);
avg_cns=avg_cns+gmsub_norm;
avg_rs = avg_rs + raw_rad;
avg_ds = avg_ds + dimension;
avg_caps_theory = avg_caps_theory + 1/capacity_theory;
end

avg_cns = avg_cns/num_to_get_geom;
avg_rs = avg_rs/num_to_get_geom;
avg_ds = avg_ds/num_to_get_geom;
avg_caps_theory = 1/((avg_caps_theory)/num_to_get_geom);
capacity_theory = 1/width;

cell_spk_counts = {};

for t = 1:length(thetas)
  cell_spk_counts{t} = counts((t-1)*num_trials+1:t*num_trials,:);
end
options = struct();
options.seed0 = 102194;
options.flag_NbyM = 0;
options.n_rep = 1;

capacity_sim = 0;
end