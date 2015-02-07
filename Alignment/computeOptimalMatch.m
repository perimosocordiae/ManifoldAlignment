%Consider the best match of two distance matrices
function dist=computeOptimalMatch(matrix1, matrix2)
  k=size(matrix1,1);
  P=getAllCombinations(k-1);
  dist=10000000000;
  for i=1:size(P,1)
      tpmatrix2=switchRowandCollumn(matrix2,P(i,:));
      tpdist=compareDistanceMatrices(matrix1, tpmatrix2);
      if tpdist<dist
          dist=tpdist;
      end
  end
end



%generate a new  matrix from an old one, the column
%and row of the new one are from the old one.
function matrix2=switchRowandCollumn(matrix1, order)
%order n*1 matrix
%matrix1 n*n matrix
  n=size(matrix1,1);
  for i=1:n
      matrix2(i,:)=matrix1(order(i),:);
  end
  matrix1=matrix2;
  for i=1:n
      matrix2(:,i)=matrix1(:,order(i));
  end
end
 


%Given k, list all the combinations of (2, 3,...,k+1)
%the first element is always 1.
function P=getAllCombinations(k)
  %List all permutatoins of v
  v=size(k);
  for i=1:k
      v(i)=i+1;
  end
  P=perms(v);
  
  %Add one stable column to P
  One=ones(size(P,1),1);
  P=[One P];
end



%Compute the smallest distance between two distance matrices.
function distance= compareDistanceMatrices(X, Y)
  %change Y;
  k=trace(Y'*X)/trace(Y'*Y);
  dist1=norm(X-k*Y);
  %change X;
  k=trace(X'*Y)/trace(X'*X);
  dist2=norm(Y-k*X);
  %the result is the smallest of the two above.
  distance=min(dist1, dist2);
end