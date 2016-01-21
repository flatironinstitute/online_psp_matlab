function err=compute_reconstruction_error(V_exact,V_approx)
%{ 
Compute quickly and efficiently the reconstruction errors given a subspace approximation. 

the input vector are given as d x q where d is the dimensionality and q is
the number of principal components

Parameter
---------
V_exact,V_approx

Return
------
err: the error ||VV_exact*V__exact' - V_approx*V_approx'||/||VV_exact*V__exact'||

%}

F=V_approx;
Vb=V_exact;

q=size(V_exact,2);

A=Vb'*F;
B=F'*F;

err=(q-2*trace(A*A')+trace(B^2))/q;
