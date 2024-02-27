function [c_opt,R_opt,t_opt] = NonMinimalSolver_CR(Scene,Shapes,K,lambda,best_set)

best_size=numel(best_set);

if best_size<6
    
disp('IMOT fails: Not enough correspondence!');
c_opt=zeros(K,1);
R_opt=eye(3);
t_opt=zeros(3,1);
return
    
end

Scene_mean=mean(Scene(:,best_set),2);

Shapes_mean=zeros(3,K);
for i=1:K
    Shapes_mean(:,i)=mean(Shapes(:,best_set,i),2);
end


Scene_=zeros(3*best_size,1);Shapes_=zeros(3*best_size,K);
for i=1:best_size
    Scene_(3*i-2:3*i)=Scene(:,best_set(i))-Scene_mean;
    Med=Shapes(:,best_set(i),:);
    Shapes_(3*i-2:3*i,:)=reshape(Med(:),3,[])-Shapes_mean;
end

B_=Shapes_;
y_=Scene_;
Y=reshape(y_,3,[]);
H_=2*(B_'*B_+lambda*eye(K));
One=ones(K,1);
G=inv(H_)-(inv(H_)*One*One'*inv(H_))/((One'*inv(H_)*One));
g=(inv(H_)*One)/(One'*inv(H_)*One);


M=[2*B_*G*B_'-eye(3*best_size);2*sqrt(lambda)*G*B_'];

h=[B_*g;g];

P=zeros(9,9);
P(1,1)=1;P(2,4)=1;P(3,7)=1;
P(4,2)=1;P(5,5)=1;P(6,8)=1;
P(7,3)=1;P(8,6)=1;P(9,9)=1;

Q=[h'*h,h'*M*kron(Y',eye(3))*P;
   (h'*M*kron(Y',eye(3))*P)',P'*kron(Y,eye(3))*M'*M*kron(Y',eye(3))*P];

mpol q 4;
R = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2) - q(3)*q(4)), 2*(q(1)*q(3) + q(2)*q(4));
        2*(q(1)*q(2) + q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(2)*q(3)-q(1)*q(4));
        2*(q(1)*q(3) - q(2)*q(4)), 2*(q(2)*q(3)+q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^2];
r=vec(R);

f_=[1;r]'*Q*[1;r];

P= msdp(min(f_),...
   q'*q==1);
[~,~] = msol(P);

q_mat  = double(mom(q*q'));
[V,~] = eig(q_mat);
R_opt = quaternion2rot(V(:,4));

c_opt=2*G*B_'*kron(eye(best_size),R_opt')*y_+g;

t_opt=Scene_mean-R_opt*sum(c_opt'.*Shapes_mean,2);   


end

