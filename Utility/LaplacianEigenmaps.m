%Laplacian eigenmaps.
function [g1]=  LaplacianEigenmaps(N1, epsilon)
%N1: n1*k matrix. k nearest neighbours for each example.
    
    %~~~Default Parameters~~~
    m=2000;    %max dimensionality of the new space.

    
    %construct weight matrix for each domain
    n1=size(N1,1);
    K=size(N1,2);

    
    W1=sparse(n1,n1);
    for i=1:n1; 
        for j=1:K;
            W1(i,N1(i,j))=1;
        end
    end
    W1=0.5*(W1+W1');
    
    W1=sparse(W1);
    
    D=sum(W1);
    D=sqrt(D);
    N=size(W1,1);
    D21=sparse(N,N);
    I=sparse(N,N);
    for i=1:max(size(D,1), size(D,2));
        I(i,i)=1;
    	if D(i)==0
    		D21(i,i)=1;
        else
            D21(i,i)=1/D(i);
    	end
   end
   
   W1=sparse((D21))*W1*sparse((D21));

   %~~~eigen decomposition~~~
   [ev, ea]=eig(full(I-W1));


    %sorting ea by ascending order
    ea=diag(ea);
    [x, index]  =sort(ea);
    ea =ea(index); ev=ev(:,index);
    for i=1:size(ev,2)
        ev(:,i)=ev(:,i)/norm(ev(:,i));
    end

    %some eigenvalues might be close to 0, and should be filted out
    for i=1:size(ea);
        if ea(i)>epsilon 
            break;
        end
    end
    start=i;

    %~~~compute mappings~~~
    if m>size(ev,2)-start+1 m=size(ev,2)-start+1; end
    g1=ev(1:n1,start:m+start-1);
end


