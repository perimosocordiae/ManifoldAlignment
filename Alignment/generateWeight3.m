function W= generateWeight3(Input1, Input2, k)
%
%Input1: a P*M matrix representing manifold 1, where P=#Features and M=#Examples;
%Input2: a Q*N matrix representing manifold 2, where Q=#Features and N=#Examples;
%W: M*N Matrix, where W(i,j)= similarity of Input1(i) and Input2(j)
%k: k in kNN

%get the size parameters of the input
P=size(Input1,1);
M=size(Input1,2);
Q=size(Input2,1);
N=size(Input2,2);

%decompose the input structures into substructures
%decompose:  regular 1D distance
%decompose2: regular 1D similarity
%decompose3: distance matrix
%decompose4: similarity matrix
[Dist1, Pos1]=decompose3(Input1, k); %each column is a substructure
[Dist2, Pos2]=decompose3(Input2, k); %each column is a substructure

%generate the weight matrix
%scaler=sum(sum(SI1))/(M*M);
W=zeros(M,N);
for i=1:M
    for j=i:N
        W(i,j)=computeOptimalMatch(squeeze(Dist1(i,:,:)), squeeze(Dist2(j,:,:)));
        W(j,i)=W(i,j);
    end
    i
end
maxmaxw=max(max(W));
W=maxmaxw-W;

%only keep the k nearest neighbors
for i=1:M
    [sorted,index]=sort(W(i,:),'descend');
    W(i,:)=W(i,:).*(W(i,:)>=sorted(k));
end
W=(W+W')*0.5;

end
