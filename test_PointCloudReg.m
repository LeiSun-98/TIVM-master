clc;
clear all;
close all;

% add dependencies (solvers)
addpath lib;

% Load data
load('data\data_PointCloudReg.mat');

% Select mode: 
% 1 for known inlier-noise statistics
% 0 for unknown inlier-noise statistics
mode_flag=1;

% Set inlierthreshold
inlier_thres=5*noise;

% Whether to plot the correspondences (1 for yes, 0 for no)
if_plot=1;

for itr_outlier=1:9

outlier_ratio=outliers(itr_outlier+1);

R_gt=cell2mat(store_R(itr_outlier+1));
t_gt=cell2mat(store_t(itr_outlier+1));
pts_3d=cell2mat(store_pts_3d(itr_outlier+1));
pts_3d_=cell2mat(store_pts_3d_(itr_outlier+1));

disp(['Noise Level: ',num2str(noise)]);
disp(['Correspondence Number: ',num2str(n_ele)]);
disp(['Outlier Ratio: ',num2str(outlier_ratio*100),'%']);

% Plot correspondences
if if_plot
    plot_correspondence(pts_3d,pts_3d_,outlier_ratio);
end

%% TBVM

tic;

re=zeros(1,n_ele);thres=zeros(1,100);best_set=1:n_ele;
ostu_itr=2;max_itr_TBVM=100;mean_re_last=0;check_converge=0;

for NM_itr=1:max_itr_TBVM

[R_opt,t_opt] = NonMinimalSolver_PCR(pts_3d,pts_3d_,best_set);

for i=1:n_ele
    re(i)=norm(R_opt*pts_3d(i,:)'+t_opt-pts_3d_(i,:)');
end

Upper=max(re);

unit=Upper/300;

unit_num=floor(Upper/unit)+1;

votes=zeros(1,unit_num);

for i=1:n_ele
    x=floor(re(i)/unit)+1;
    votes(x)=votes(x)+1;
end

coe=zeros(1,1);
for ree=1:ostu_itr

serial=1:numel(votes);
p = votes' / sum(votes);
omega = cumsum(p);
mu = cumsum(p .* (serial)');
mu_t = mu(end);

Gvariance=sum((serial'-mu_t).^2.*p);

sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

bcvariance=sigma_b_squared/Gvariance;

[max_bc,vote_idx]=max(bcvariance);

votes=votes(1:vote_idx);

coe(ree)=max_bc;

end

thres(NM_itr)=unit*vote_idx;

best_set_last=best_set;
co=0;
best_set=ones(1,1);
for i=1:n_ele
    if re(i)<=thres(NM_itr)
        co=co+1;
        best_set(co)=i;
    end
end

if (check_converge>0 && abs(mean(re)-mean_re_last)<=1e-3*mean_re_last) 
    break
else
    check_converge=0;
end

if mode_flag==1 && thres(NM_itr)<=2*inlier_thres
    break
end

if NM_itr>1 && abs(thres(NM_itr)-thres(NM_itr-1))<unit  
    ostu_itr=ostu_itr+1;
    check_converge=1;
    mean_re_last=mean(re);
end


end

itr_additional=0; 

if mode_flag==1

    best_set=ones(1,1);count=0;
    for i=1:n_ele
        if re(i)<=inlier_thres
            count=count+1;
            best_set(count)=i;
        end
    end

    [R_opt,t_opt] = NonMinimalSolver_PCR(pts_3d,pts_3d_,best_set);

    itr_additional=itr_additional+1;

end

disp(['Rotation Error [deg]: ',num2str(compAngErr(R_opt,R_gt)*180/pi)]);
disp(['Translation Error [m]: ',num2str(norm(t_opt-t_gt'))]);
disp(['Number of Iterations:',num2str(NM_itr+itr_additional)]);
disp(' ');

end

