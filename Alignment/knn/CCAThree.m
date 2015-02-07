%CCA: three domains
function [cca1, cca2, cca3]=  CCAThree(X1, X2, X3, N1, N2, N3, W12, W13, W23, epsilon, mu)
%X1: P1*M1 matrix, M1 examples in a P1 dimensional space.
%X2: P2*M2 matrix
%X3: P3*M3 matrix
%N1: M1*k matrix. k nearest neighbours for each example.
%N2: M2*k matrix
%N3: M3*k matrix
%W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
%W13: M1*M3 sparse matrix
%W23: M2*M3 sparse matrix
%epsilon: precision. 
%\mu: NOT USED
   
    %~~~Default Parameters~~~
    m=2000;    %max dimensionality of the new space.
    
    n1=size(N1,1); 
    n2=size(N2,1); 
    n3=size(N3,1);
    K=min (size(N1,2), size(N2,2));
    K=min (K,size(N3,2));
    W1=sparse(n1,n1); 
    W2=sparse(n2,n2); 
    W3=sparse(n3,n3);
    
    
    W=[0*W1 W12 W13; W12' 0*W2 W23; W13' W23' 0*W3]; 
    W=sparse(W);
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

    %~~~size of X and Y~~~
    P1=size(X1,1); M1=size(X1,2);  
    P2=size(X2,1); M2=size(X2,2);  
    P3=size(X3,1); M3=size(X3,2);  


    %create W, D, L, Z
    Z=[X1 zeros(P1,M2+M3); zeros(P2, M1) X2 zeros(P2, M3); zeros(P3, M1+M2) X3];
    Z=sparse(Z);

    %Create T, Tplus (T^+)
    [u, s, v]=svd(full(Z*Z'));
    F=u*sqrt(s);
    Fplus=pinv(F);
    
    clear u s v;
    save data.mat;
    
    T=Fplus*Z*(I-W)*Z'*Fplus'; 

 
    %~~~eigen decomposition~~~
    T=0.5*(T+T');
    [ev, ea]=eig(full(T));
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
    if m>size(ev,2)-start+1 m=size(ev,2)-start+1; end
    cca1=ev(1:P1,start:m+start-1);
    cca2=ev(P1+1:P1+P2, start:m+start-1);
    cca3=ev(P1+P2+1:P1+P2+P3, start:m+start-1);
end


