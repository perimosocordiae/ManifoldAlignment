function Wx=createAllConnectedGraph(X, delta)

    maxx=max(max(abs(X))); 
    Wx=L2_distance(X/maxx, X/maxx, 0); 
    Wx=exp(-Wx.^2/delta^2);
    Dxsqrt=1./sqrt(sum(Wx)); 
    Dxsqrt=diag(Dxsqrt); 
    Wx=Dxsqrt*Wx*Dxsqrt; 
end    
    