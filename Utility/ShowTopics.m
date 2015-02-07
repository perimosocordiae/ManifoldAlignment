function [str, vtr]=ShowTopics (cols, absolute, wl2)

Topics=full(cols);
m=size(Topics,2);

% Print out topics to str
for i=1:m
  if absolute==1
     [value, id]=sort(abs(Topics(:,i)), 'descend');
  else 
     [value, id]=sort(Topics(:,i), 'descend');
  end

  for j=1:10;
    str(i,j)  =wl2(id(j));
    vtr(i,j)=value(j)*value(j);
  end
end

end