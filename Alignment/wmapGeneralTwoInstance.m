%Instance-Level Manifold Projections. Two domains
function [g1, g2]=  wmapGeneralTwoInstance(X1, X2, W1, W2, W12, epsilon, mu)
%X1: NOT USED. %P1*M1 matrix, M1 examples in a P1 dimensional space.
%X2: NOT USED. %P2*M2 matrix
%N1: M1*M1 matrix. weight matrix for each domain.
%N2: M2*M2 matrix
%W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
%epsilon: precision. 
%\mu: used to balance two goals: matching corresponding pairs and preserving manifold topology. 
    
    %~~~Default Parameters~~~
    global DimLatentSpace;    %max dimensionality of the new space.
    
    %construct weight matrix for each domain
    n1=size(W1,1); 
    n2=size(W2,1);
        
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
     
%    [ev, ea]=eig(full(I-W));

    [ev, ea]=eigs(full(I-W),m,'SM');


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
    if DimLatentSpace>size(ev,2)-start+1 DimLatentSpace=size(ev,2)-start+1; end
    g1=ev(1:n1,start:m+start-1);
    g2=ev(n1+1:n1+n2, start:m+start-1);
end


