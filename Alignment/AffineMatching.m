function P = AffineMatching (Xtrain, Ytrain);

%Xtrain: p*m matrix: m examples, each of which is represented by p features.
%Ytrain: p*m matrix

 
%Here, we want to find P such that X=PY, so P=X/Y;
P=Xtrain/Ytrain;

return;