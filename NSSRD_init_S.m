function  S=NSSRD_init_S(fea,m)

[N,dim]=size(fea);
options.PCARatio=0.75;%0.75
[eigvector, eigvalue] = PCA(fea,options);
fea=fea*eigvector;
label=litekmeans(fea,m,'MaxIter',100,'Replicates',10);
post=zeros(N,m);
for i=1:N
    %  post(i,label_new(i))=1;%将Pxi矩阵转化为0-1矩阵
    post(i,label(i))=1;
end
S=post';