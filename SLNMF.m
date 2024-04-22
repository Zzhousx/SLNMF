function [X_new]=SLNMF(class_num,X,W,H,G,alpha,beta,beta1,gamma,r,l,NIter)
[n,d]=size(X);
U=eye(d);
P=eye(n);
Cn=eye(n)-ones(n)/n;
F=initfcm(n,r);
M = rand(d, class_num);
for iter=1:NIter
    for j=1:class_num
        aa = 0;
       
        for i=1:n

            aa = aa + F(i,j)*X(i,:)';
      
        end
        M(:,j)= aa/sum(F(:,j));
    end
    for i = 1:n
        for j = 1:class_num
            distance(i,j) = mydist(X(i,:)',M(:,j));
            distance(i,j) = distance(i,j)^2;
        end
        
        ad = 0.5*(alpha*(Cn(i,:)*X*W*G) -(distance(i,:)/(2*beta)));
        
        F(i,:) = EProjSimplex_new(ad);
    end
    W=W.*((2*X'*P*X*H'+2*alpha*X'*Cn*F*G'+2*gamma*W)./(2*X'*P*X*W*H*H'+2*alpha*X'*Cn*X*W*G*G'+2*beta1*U*W+2*gamma*W*W'*W+eps));
    H=H.*((W'*X'*P*X)./(W'*X'*P*X*W*H+eps));
    G=G.*((alpha*W'*X'*Cn*F)./(alpha*W'*X'*Cn*X*W*G+eps));
    E = X-X*W*H;
    Ei = sqrt(sum(E.*E,2)+eps);%X-XWH=E的21范数
    p = 0.5./Ei;
    P = diag(p);
    Wj = sqrt(sum(W.*W,2)+eps);%W的21范数
    u = 0.5./Wj;
    U = diag(u);
end
score= sqrt(sum(W.*W,2));
[~, idx] = sort(score,'descend');
X_new = X (:,idx(1:l));
end





