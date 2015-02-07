%Compare embedding results
function result=cmpEmbedding(X, Y, f, g)
%X: p*m matrix
%Y: q*m matrix
%X(:,i)<->Y(:,i)
%f: p*d matrix. Mapping function.
%g: q*d matrix. Mapping function.


Xn=f'*X;
Yn=g'*Y;

total=0;
n=min (size(X,2), size(Y,2));
result=zeros(n,1);
x=zeros(n,1);
y=zeros(n,1);

for i=1:n;
    x=Xn(:,i);
    for j=1:n;
        y(j)=norm(x- Yn(:,j));
    end
    [v, id]=sort(y);
    
    %check where the true match is
    for k=1:n;
        if id(k)==i
            result(i)=k;
            if k==1
                total=total+1;
            end
  %          fprintf(1, '%d: %d  %d\n', i, k, total);
            break;
        end
    end
    
end
save result.mat result;

end
    
