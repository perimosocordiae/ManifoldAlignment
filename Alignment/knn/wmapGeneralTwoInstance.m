%Instance-Level Manifold Projections. Two domains
function [g1, g2]=  wmapGeneralTwoInstance(X1, X2, N1, N2, W12, epsilon, mu)
%X1: NOT USED. %P1*M1 matrix, M1 examples in a P1 dimensional space.
%X2: NOT USED. %P2*M2 matrix
%N1: M1*k matrix. k nearest neighbours for each example.
%N2: M2*k matrix
%W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
%epsilon: precision. 
%\mu: used to balance two goals: matching corresponding pairs and preserving manifold topology. 
    
    %~~~Default Parameters~~~
    m=2000;    %max dimensionality of the new space.
    
    %construct weight matrix for each domain
    n1=size(N1,1); 
    n2=size(N2,1);
    K=min (size(N1,2), size(N2,2));
    
    W1=sparse(n1,n1); 
    W2=sparse(n2,n2); 
    for i=1:n1; 
        for j=1:K;
            W1(i,N1(i,j))=1;
        end
    end
    for i=1:n2; 
        for j=1:K;
            W2(i,N2(i,j))=1;
        end
    end
  
    W1=0.5*(W1+W1');
    W2=0.5*(W2+W2');
    
    sum1=sum(sum(W1))+sum(sum(W2));
    sum2=2*(sum(sum(W12)));
    mu=mu*sum1/sum2;
    
    W=[W1 mu*W12; mu*W12' W2]; 
    W=sparse(W);
    clear W1 W2 W12;
    
    D=sum(W);
    D=sqrt(D);
    N=size(W,1);
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
   
   W=sparse((D21))*W*sparse((D21));

    save datainstance.mat;
    
     %~~~eigen decomposition~~~
    [ev, ea]=eig(full(I-W));


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
    g2=ev(n1+1:n1+n2, start:m+start-1);
end


