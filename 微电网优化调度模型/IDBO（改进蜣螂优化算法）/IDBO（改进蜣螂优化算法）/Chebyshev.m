function X=Chebyshev(N,dim,ub,lb)
k = 4;
X = rand(N,dim); 
for i = 1 : dim
    for j = 1 : N-1
        X(j+1,i) = cos(k.*acos(X(j,i)));
    end
end
for i=1:N
    X(i,:)=X(i,:).*(ub-lb)+lb;
end
end