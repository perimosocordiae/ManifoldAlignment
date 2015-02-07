function [Dist, Position]=decompose3(Input, k)
%
%Input: a P*M matrix representing manifold 1, where P=#Features and
%M=#Examples;
%k: k in kNN
%Subset: a k*M matrix representing the distance from each example to its k
%nearest neighbors.
%
M=size(Input,2);
D=L2_distance(Input, Input,0); %M*M matrix
mmax=max(max(D));

Position=zeros(M,k+1);
Dist=zeros(M,k+1,k+1);
for i=1:M
    D(i,i)=mmax;
    [value, id]=sort(D(:,i));
    Position(i,1)=i;
    for j=2:k+1
        Position(i,j)=id(j-1);
    end
    for m=1:k+1;
        for n=1:k+1;
            Dist(i,m,n)=D(Position(i,m),Position(i,n));
        end
        Dist(i,m,m)=0;
    end
end

end
    
