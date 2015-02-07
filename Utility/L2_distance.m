function d = L2_distance(a, b, df)
%   a and b are d*n matrices. n examples, each example is in R^d space.

  if (size(a,1) == 1)
    a = [a; zeros(1,size(a,2))]; 
    b = [b; zeros(1,size(b,2))]; 
  end
  
  aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
  d = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
  d = real(d); 
  if (df==1)
    d = d.*(1-eye(size(d)));
  end
