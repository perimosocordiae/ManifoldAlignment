function W=createKnnGraph(N1)

    n1=size(N1, 1);
     K=size(N1, 2);
     
    W=sparse(n1,n1); 
    
    for i=1:n1; 
        for j=1:K;
            W(i,N1(i,j))=1;
        end
    end
    W=0.5*(W+W');
end