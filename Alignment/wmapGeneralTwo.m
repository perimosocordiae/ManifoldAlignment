%Feature-level Manifold Projections. Two domains.
function [map1, map2]=  wmapGeneralTwo(X1, X2, W1, W2, W12, epsilon, mu)
%X1: P1*M1 matrix, M1 examples in a P1 dimensional space.
%X2: P2*M2 matrix
%W1: M1*M1 matrix. weight matrix for each domain.
%W2: M2*M2 matrix
%W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
%epsilon: precision. 
%\mu: used to balance two goals: matching corresponding pairs and preserving manifold topology.    
 
 
    %~~~Default Parameters~~~
    global DimLatentSpace;    %max dimensionality of the new space.
       
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

    %~~~size of X and Y~~~
    P1=size(X1,1); M1=size(X1,2);  
    P2=size(X2,1); M2=size(X2,2);  


    %create W, D, L, Z
    Z=[X1 zeros(P1,M2); zeros(P2, M1) X2];
    Z=sparse(Z);

    %Create T, Tplus (T^+)
    
%    [u, s, v]=svd(full(Z*Z'));

    [u, s, v]=svds(full(Z*Z'),DimLatentSpace,'SM');
    
    F=u*sqrt(s);
    Fplus=pinv(F);
    
    clear u s v;
    save data.mat;
    
    T=Fplus*Z*(I-W)*Z'*Fplus'; 

 
    %~~~eigen decomposition~~~
    T=0.5*(T+T');

%    [ev, ea]=eig(full(T));

    [ev, ea]=eigs(full(T),DimLatentSpace,'SM');

    clear T Z F;

    %sorting ea by ascending order
    ea=diag(ea);
    [x, index]  =sort(ea);
    ea =ea(index); ev=ev(:,index);
    ev =Fplus'*ev;
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
    map1=ev(1:P1,start:m+start-1);
    map2=ev(P1+1:P1+P2, start:m+start-1);
end


