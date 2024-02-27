clc;
clear all;
close all;

% add dependencies (solvers)
addpath lib;

% Load data
load('data\data_RotAvg.mat');

% Select mode: 
% 1 for known inlier-noise statistics
% 0 for unknown inlier-noise statistics
mode_flag=1;

% Set inlierthreshold
inlier_thres=3*noise;

for itr_outlier=1:7

outlier_ratio=outlier_ratios(itr_outlier);

R_Samples=store_R_samp{itr_outlier};

R_gt=store_R{itr_outlier};

R_input=R_Samples;

disp(['Noise Level [deg]: ',num2str(noise)]);
disp(['Outlier Ratio: ',num2str(outlier_ratio*100),'%']);


%% TIVM

tic;

n_ele=n_inliers+n_outliers;re=zeros(1,n_ele);thres=zeros(1,100);best_set=1:n_ele;
ostu_itr=2;max_itr_TBVM=100;
mean_re_last=0;check_converge=0;

for NM_itr=1:max_itr_TBVM

R_samples=R_input(best_set);
n_samples = length(R_samples);

vectors_total = zeros(9,n_samples);
for i = 1:n_samples
    vectors_total(:,i)= R_samples{i}(:);
end
s = median(vectors_total,2);
    
[U,~,V] = svd(reshape(s, [3 3]));
R = U*V.';
if (det(R) < 0)
    V(:,3) = -V(:,3);
    R = U*V.';
end

R_opt=R;

for i=1:n_ele
    re(i)=AngularError(R,R_input{i})*180/pi;
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

thres(NM_itr)=unit*vote_idx-unit;

best_set=zeros(1,1);co=0;
for i=1:n_ele
    if re(i)<=thres(NM_itr)
        co=co+1;
        best_set(co)=i;
    end
end

if check_converge>0 && abs(mean_re_last-mean(re))<=1e-3*mean_re_last 
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

    best_set=ones(1,1);co=0;
    for i=1:n_ele
        if re(i)<=inlier_thres
            co=co+1;
            best_set(co)=i;
        end
    end
    
    R_samples=R_input(best_set);

    n_samples = length(R_samples);

    vectors_total = zeros(9,n_samples);
    for i = 1:n_samples
        vectors_total(:,i)= R_samples{i}(:);
    end
    s = median(vectors_total,2);

    [U,~,V] = svd(reshape(s, [3 3]));
    R_opt = U*V.';
    if (det(R_opt) < 0)
        V(:,3) = -V(:,3);
        R_opt = U*V.';
    end

    itr_additional=1;

end

disp(['Rotation Error [deg]: ',num2str(compAngErr(R_opt,R_gt)*180/pi)]);
disp(['Number of Iterations:',num2str(NM_itr+itr_additional)]);
disp(' ');

end


