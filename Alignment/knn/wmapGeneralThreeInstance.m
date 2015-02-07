%Instance-Level Manifold Projections. Three domains.
function [g1, g2, g3]=  wmapGeneralThreeInstance(X1, X2, X3, N1, N2, N3, W12, W13, W23, epsilon, mu)
%X1: NOT USED. P1*M1 matrix, M1 examples in a P1 dimensional space.
%X2: NOT USED. P2*M2 matrix
%X3: NOT USED. P3*M3 matrix
%N1: M1*k matrix. k nearest neighbours for each example.
%N2: M2*k matrix
%N3: M3*k matrix
%W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
%W13: M1*M3 sparse matrix
%W23: M2*M3 sparse matrix
%epsilon: precision. 
%\mu: used to balance two goals: matching corresponding pairs and preserving manifold topology.   
    
    %~~~Default Parameters~~~
    m=2000;    %max dimensionality of the new space.
    
    %construct weight matrix for each domain
    n1=size(N1,1); 
    n2=size(N2,1); 
    n3=size(N3,1);
    K=min (size(N1,2), size(N2,2));
    K=min (K,size(N3,2));
    W1=sparse(n1,n1); 
    W2=sparse(n2,n2); 
    W3=sparse(n3,n3);

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
    for i=1:n3; 
        for j=1:K;
            W3(i,N3(i,j))=1;
        end
    end
    W1=0.5*(W1+W1');
    W2=0.5*(W2+W2');
    W3=0.5*(W3+W3');
    
    sum1=sum(sum(W1))+sum(sum(W2))+sum(sum(W3));
    sum2=2*(sum(sum(W12))+sum(sum(W13))+sum(sum(W23)));
    mu=mu*sum1/sum2;
    
    W=[W1 mu*W12 mu*W13; mu*W12' W2 mu*W23; mu*W13' mu*W23' W3]; W=sparse(W);
    clear W1 W2 W3 W12 W13 W23;
    
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
    g3=ev(n1+n2+1:n1+n2+n3, start:m+start-1);
end


