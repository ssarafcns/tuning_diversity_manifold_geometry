function [raw_rad, dimension, gmsub_norm, width, capacity_theory, capacity_sim, avg_cns, avg_rs, avg_ds, avg_caps_theory, phis_angles,lambda, mean_dots_noise_cent,exp_rest_avg, theta_rests,total_num_supp_vec, corrs_per_dist] = manifold_heterogeneity_FR_parallel_preexist_2(data_mat,pa_given)

num_trials = size(data_mat,2);
num_thetas = size(data_mat,3);
N = size(data_mat,1);
counts = zeros(num_thetas*num_trials, N);
sub_counts = zeros(length(pa_given)*num_trials,N);
for i = 0:num_thetas-1
counts(i*num_trials+1 : (i+1)*num_trials,1:N) = data_mat(:,:,i+1)';

if ~isempty(find(pa_given==i+1))
    ind_in_given = find(pa_given==i+1)-1;
    sub_counts(ind_in_given*num_trials+1 : (ind_in_given+1)*num_trials,1:N) = data_mat(:,:,i+1)';
end
end

centroids = zeros(num_thetas,N);
centroids_unglob_unnorm = zeros(num_thetas,N);
centroids_unnorm = zeros(num_thetas,N);
mean_dots_noise_cent_vec = zeros(1,num_thetas);
exp_rest_vec = zeros(1, num_thetas);
centroids_rest = zeros(num_thetas,N);

exp = {};
exp_rest = {};
mani_axes = {};
mani_axes_rest = {};
global_mean = mean(counts,1);
d = reshape(data_mat, [N, num_thetas*num_trials]);
d= d';
theta_rests = {};
total_num_supp_vec = zeros(num_thetas,1);

for theta_num = 0:num_thetas-1

  centroids(theta_num+1,:) = mean(counts(theta_num*num_trials+1:theta_num*num_trials+num_trials,:),1);
  centroids_unglob_unnorm(theta_num+1, :) = centroids(theta_num+1,:);
  centroids(theta_num+1,:) = centroids(theta_num+1,:)-global_mean;
  centroids_unnorm(theta_num+1,:) = centroids(theta_num+1,:);
  centroids(theta_num+1,:) = centroids(theta_num+1,:)/norm(centroids(theta_num+1,:));
  
  [coeffs, ~,~, ~, exp{theta_num+1},~] = pca(counts(theta_num*num_trials+1:theta_num*num_trials+num_trials,:));
  ma_hold = zeros(length(coeffs(1,:)),N);
  for ma=1:length(coeffs(1,:))
  ma_hold(ma,:) = ( coeffs(:,ma)/norm(coeffs(:,ma)) )';
  end
  mani_axes{theta_num+1} = ma_hold;


end
exp_rest_avg = mean(exp_rest_vec);
pa_to_use = [pa_given];
num_compute = length(pa_to_use);
C_albert_mat = zeros(72,72);
phis_angles = zeros(num_thetas,num_compute); 
for pa=pa_to_use
 pa_thetas = mod(pa+[ 1], num_thetas);
 pa_thetas(pa_thetas==0)= num_thetas;
 non_spont_neur = 1:N;
 centroid_both = [centroids(pa,:); centroids(pa_thetas(1),:)];

 inn_prod_cen = centroids*centroids';

 phis_angles(:,:) = inn_prod_cen;
end
unit_centroid_mat = zeros(length(pa_given),N);

C_albert_mat = centroids(pa_given,:)*centroids(pa_given,:)';
es = eig(C_albert_mat);
es = sort(es, 'ascend');

corrs_per_dist = 0;

lambda= [];
mean_dots_noise_cent = [];
count_pairs = 0;
num_axes = min(num_trials,N)-1;
exp_tot_mani_1 = 0;
if length(pa_given) == 5
    other_manis = [-2:2];
elseif length(pa_given) == 12
    other_manis = [-5:6];
else
    other_manis = [-17:18];
end
for theta_num = pa_to_use
    l=0;
    mean_dots_noise_c = 0;
    mean_dots_c = [];
    exp_now = exp{theta_num};
    mani_axes_now = mani_axes{theta_num};
    l_overlap = [];
  
    for other_mani = -17:18

        other_mani_mod = mod(theta_num+other_mani,num_thetas);
         if other_mani_mod == 0
             other_mani_mod = num_thetas;
         end
         exp_other = exp{other_mani_mod};
         mani_axes_other = mani_axes{other_mani_mod};
         l_hold = 0;
         unit_cent = centroids(theta_num, :)/norm(centroids(theta_num,:));

    for m_axis_1 = 1:num_axes

      exp_tot_mani_2 = 0;
      for m_axis_2 = [m_axis_1]
        l_hold = l_hold + ( - acos(abs(dot(mani_axes_now(m_axis_1,:),mani_axes_other(m_axis_2,:)))))*exp_now(m_axis_1)/sum(exp_now);
        
      end
    end
    lambda(theta_num, other_mani_mod) = l_hold;

    end

    for m_axis_1 = 1:num_axes
      mean_dots_noise_c= mean_dots_noise_c + abs(dot(mani_axes_now(m_axis_1,:), unit_cent))*exp_now(m_axis_1)/sum(exp_now);
    end
    mean_dots_noise_cent = [mean_dots_noise_cent, mean_dots_noise_c];
end




avg_cns = [];
avg_rs = [];
avg_ds = [];
avg_caps_theory = 0;
num_to_get_geom = 2;
new_pa_to_use = pa_to_use;

for t=new_pa_to_use

other_theta_c = mod(t+1, num_thetas);
if other_theta_c == 0
    other_theta_c =num_thetas;
end
mean_two_centers = mean([centroids_unglob_unnorm(t,:);centroids_unglob_unnorm(other_theta_c,:)], 1);

[~,dimension, ~,gmsub_norm,raw_rad,~,~, ~, width, capacity_theory] = calculateManifoldDimensionRadius_Gaussian(counts((t-1)*num_trials+1:t*num_trials,:),global_mean, num_to_get_geom);

[~,d2, ~,gmsub_norm_2,r2,~,~, ~, ~, ~] = calculateManifoldDimensionRadius_Gaussian(counts((other_theta_c-1)*num_trials+1:other_theta_c*num_trials,:),global_mean, num_to_get_geom);

avg_cns = [avg_cns, gmsub_norm];
avg_rs = [avg_rs, raw_rad];
avg_ds = [avg_ds, dimension];

avg_caps_theory = avg_caps_theory + 1/capacity_theory;
end

avg_caps_theory = 1/(avg_caps_theory/length(new_pa_to_use));

for t = 1:num_thetas
  cell_spk_counts{t} = counts((t-1)*num_trials+1:t*num_trials,:);
end
options = struct();
options.seed0 = 102194;
options.flag_NbyM = 0;
options.n_rep = 1;


capacity_sim = 0;
end