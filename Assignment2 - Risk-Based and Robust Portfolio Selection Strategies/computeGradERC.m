
function gval = computeGradERC (x)

global Q
  
  n = size(Q,1) ;  

  if(size(x,1)==1)
     x = x';
  end
  
  % Insert your gradiant computations here
  % You can use finite differences to check the gradient

  
  gval = zeros(n,1); 
  y = x .* (Q*x);
    
    Q_x = zeros(n,1);  % y = x.Qx, here we calculate Qx to help gradient formulation
    for i = 1:n
        Q_x(i) = Q(i,:) * x;
    end
    % so now we have a matrix of x*1  --> Qx, each of the entry will help
    % construct the partial derivatives in the gradient function  eg: Q_x(1) = a11x1 + a12x2 + a13x3
    
    % for each row of the gradient matrix, it consists of sum of the
    % products between (Yi-Yj) and (partial derivative of Ya on Xb - partial derivatiev of Yc on Xd)
    % to figure out what i,j,a,b,c,d are, we have found the pattern and construct the model using nested loop.
    
    for k = 1:n % for each row of the gradient matix
        for i = 1:n % i is the subscript coefficient of y in the first partial derivative
            for j = i+1:n % j is the subscript coefficient of y in the second partial derivative
                if i == k % when starting Y(i) corresponds to kth gradient eg: k = 1, i = 1, j = 2
                    % the first partial derivative is eg: 2a11x1 + a12x2 + a13x3 = Q_x(1) + a11x1
                    % the second partial derivative is eg: a21x2                 
                    pd1 = Q_x(i) + Q(k,k)*x(k);
                    pd2 = Q(j,i)*x(j);
                    gval(k) = gval(k)+ (y(i)-y(j))*(pd1-pd2);   
                elseif j == k % when ending Y(i) corresponds to kth gradient eg: k = 2, i = 1, j =2
                    % the first partial derivative is eg: a12x1
                    % the second partial derivative is eg: a11x1 + 2a12x2 + a13x3 = Q_x(2) + a12x2
                    pd1 = Q(i,j)*x(i);
                    pd2 = Q_x(j) + Q(k,k)*x(k);
                    gval(k) = gval(k)+ (y(i)-y(j))*(pd1-pd2);
                else % for the rest of the cases, partial derivatives are simple
                    pd1 = Q(i,k)*x(i);
                    pd2 = Q(j,k)*x(j);
                    gval(k) = gval(k)+ (y(i)-y(j))*(pd1-pd2);
                end
            end
        end
    end
  gval = 2*2*gval; %times back due to Cplex 1/2 X*QX
  
  
% forward approximation
%Step Sizes
%   h=0.00001; step size the smaller the better
%   for k = 1:n
%       x2 = x;
%       x2(k) = x2(k)+h;
%       y2 = x2 .* (Q*x2);     
%       for i = 1:n
%           for j = i+1:n
%               xij = (y2(i) - y2(j)).^2 - (y(i) - y(j)).^2;
%               gval(k) = gval(k) + xij;
%           end
%       end
%   end 
%   gval = 2*gval/h;


  
  
 
end

%% Validation with finite difference methods
% function f = obj(x)
%   f  = x .* (Q*x);
% end
%
% % finite difference methox
% h = 0.0001; % step size, the smaller the better
% for i = 1:n
%     for j = i+1: n 
%          xij = (obj(x(i)+h) - obj(x(j)+h)).^2 - (obj(x(i)) - obj(x(j)).^2;
%     end
%     gval(i) = gval(i) + xij;
% end
% gval = 2*gval/h;


