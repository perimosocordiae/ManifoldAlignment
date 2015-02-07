function [Q, k, X_mean, Y_mean] = Procrustes (X, Y);

%Procrustes alignment.

%INPUT:
%X: P*M matrix. M examples, each of which is in a R^P space.
%Y: P*M matrix. M examples, each of which is in a R^P space.

%OUTPUT:
%k     : rescaling factor;
%Q     : rotation.
%X_mean: mean of X.
%Y_mean: mean of Y.


%#examples
M=min (size(X, 2), size(Y,2));

%mean
X_mean= mean(X');
Y_mean= mean(Y');

%
for i=1:M;
  X(:,i)=X(:,i)-X_mean';
  Y(:,i)=Y(:,i)-Y_mean';
end
 
%Procrustes alignment
[u, s, v]=svd(Y*X');
Q=u*v';
k=trace(s)/(trace(Y*Y'));

return;
