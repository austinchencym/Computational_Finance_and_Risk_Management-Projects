function fval = computeObjERC (x)

global Q

  n = size(Q,1) ;  

  if(size(x,1)==1)
     x = x';
  end
  
  y = x .* (Q*x) ; 
  
  fval = 0 ; 
  
  for i = 1:n
    for j = i+1:n
      xij  = y(i) - y(j) ; 
      fval = fval + xij*xij ; 
    end 
  end
  
  fval = 2*fval ;     
  
end
