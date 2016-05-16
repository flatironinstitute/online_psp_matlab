function y = NSM_neuraldynMat(Wx,M,tol,init)


y = init;

error = 1;
prev_error = 1;
k = length(M(1,:));

count = 0;

while error > tol
    
   y_old = y;
    
   for i = 1:k
       
      dum = Wx(i)-M(i,:)*y;
      y(i) = max(dum,0);
   end
   
   error = 0;
   dum = abs(y_old - y)./abs(y_old+0.001);
   
   error = max([error; dum]);
   prev_error = error;
   
   if error == 0
        if prev_error == 0
            error = 0;
        else 
            error = 1;
        end
       
   end
   
   count = count+1;
   
   if count == 2000
       disp(['******* ********* did not converge: err,' num2str(error)])
       error = 0;
       
   end
end

%count

end