%LPP
function [map1]=  LPP(X1, N1, epsilon)
%X1: P1*n1 matrix, M1 examples in a P1 dimensional space.
%N1: n1*k matrix. k nearest neighbours for each example.

    %~~~Default Parameters~~~
    
%    m=200;    %max dimensionality of the new space.
global  DimLatentSpace; 
global SVDResolution; 

    K=size(N1,2);
    
    %construct weight matrix for each domain
    n1=size(N1,1);
    
    W1=sparse(n1,n1);
    for i=1:n1; 
        for j=1:K;
            W1(i,N1(i,j))=1;
        end
    end
    W1=0.5*(W1+W1');
    
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

    %~~~size of X and Y~~~
    P1=size(X1,1); 

    %Create T, Tplus (T^+)
    %[u, s, v]=svd(full(X1*X1'));
    
    
%    fprintf('Finding %d singular vectors of X1*X1^T matrix of size %d x
%    %d\n', DimLatentSpace,size(X1,1),size(X1,1));

    fprintf('Finding singular vectors of X1*X1^T matrix of size %d x %d\n',size(X1,1),size(X1,1));
        
        
%    [u, s, v]=svds(full(X1*X1'),DimLatentSpace); 

% use faster qsvd routine -- 30 day trial! 
   [u, s, v]=qsvd(full(X1*X1'),SVDResolution); 
   
   fprintf('Found %d singular vectors upto resolution %f\n', size(u,2), SVDResolution); 
   
   
    F=u*sqrt(s);
    Fplus=pinv(F);
    
    clear u s v;
    
    T=Fplus*X1*(I-W1)*X1'*Fplus'; 
  
    %~~~eigen decomposition~~~
    T=0.5*(T+T');
   % [ev, ea]=eigs(full(T),DimLatentSpace,'SM');
    % [ev, ea]=eigs(full(T),DimLatentSpace,'SM');
    
     % [ev, ea]=eig(full(T)); % since T is of dim DimLatentSpace
     
    fprintf('Finding eigenvectors of T matrix of size %d x %d\n',size(T,1),size(T,2)); 
     
     if size(T)==min(DimLatentSpace,size(T,1))
        [ev, ea]=eig(full(T)); % since T is of dim DimLatentSpace
     else  [ev, ea]=eigs(T,DimLatentSpace,'SM');
     end;

   
    clear T F;
% 
%     %sorting ea by ascending order
    ea=diag(ea);
    [x, index]  =sort(ea);
    ea =ea(index); ev=ev(:,index);
    ev =Fplus'*ev;
    for i=1:size(ev,2)
        ev(:,i)=ev(:,i)/norm(ev(:,i));
    end

%     %some eigenvalues might be close to 0, and should be filted out
    for i=1:size(ea);
        if ea(i)>epsilon 
            break;
        end
    end
    start=i;
% 
%     %~~~compute mappings~~~
    if DimLatentSpace>size(ev,2)-start+1 DimLatentSpace=size(ev,2)-start+1; end
    map1=ev(1:P1,start:DimLatentSpace+start-1);
end


