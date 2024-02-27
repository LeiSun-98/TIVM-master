function [R_opt,t_opt] = NonMinimalSolver_PCR(pts_3d,pts_3d_,best_set)

q_=zeros(3,1);
p_=zeros(3,1);

best_size=length(best_set);

if best_size<=1
    R_opt=eye(3);
    t_opt=zeros(3,1);
    return
end

for i=1:best_size

q_(1)=q_(1)+pts_3d_(best_set(i),1);
q_(2)=q_(2)+pts_3d_(best_set(i),2);
q_(3)=q_(3)+pts_3d_(best_set(i),3);

p_(1)=p_(1)+pts_3d(best_set(i),1);
p_(2)=p_(2)+pts_3d(best_set(i),2);
p_(3)=p_(3)+pts_3d(best_set(i),3);

end

p_=p_/best_size;
q_=q_/best_size;

s_best=1;

H=zeros(3,3);
for i=1:best_size
    H=H+(pts_3d(best_set(i),:)'-p_)*(pts_3d_(best_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-s_best*R_opt*p_;


end

