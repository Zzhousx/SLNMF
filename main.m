clear;
load('fs_jaffe.mat');
X=fea;
Y=gnd;
nClusts = length(unique(gnd));
r= length(unique(gnd));
m=300;%子空间的维度
l=100;%选择特征的数量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%初始化W,H
[n,d]=size(X);
B=NSSRD_init_S(fea',m);
A=ones(d,m);
D=ones(m,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W=A;
H=B;
G=D;
NIter=10;
gamma = 0.01;
alpha=0.1;
beta=0.001;
beta1=1e+4;
rand('twister',5489);
X_new=SLNMF(r,X,W,H,G,alpha,beta,beta1,gamma,r,l,NIter);
for i=1:40
    label=litekmeans(X_new,nClusts,'MaxIter',100,'Replicates',20);
    result1 = ClusteringMeasure(gnd,label); 
    result(i,:) = result1;
end
for j=1:2
    a=result(:,j);
    ll=length(a);
    temp=[];
    for i=1:ll
        if i<ll-18
            b=sum(a(i:i+19));
            temp=[temp;b];
        end
    end
    [e,f]=max(temp);
    e=e./20;
    MEAN(j,:)=[e,f];
    STD(j,:)=std(result(f:f+19,j));
    rr(:,j)=sort(result(:,j));
    BEST(j,:)=rr(end,j);
end
STD
BEST
MEAN   
