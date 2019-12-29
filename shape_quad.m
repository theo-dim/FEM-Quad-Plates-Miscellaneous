%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shape_Quad.m - 1/12/16                                   %
% author: Tehila Stone | Theo Dimitrasopoulos              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N, B, detJ] = shape_quad(r, s, xn, ien, nen) 

x = zeros(1,nen);
y = zeros(1,nen);
    for i= 1:nen
        x(1,i) = xn(1, ien(i)); 
        y(1,i) = xn(2, ien(i)); 
    end
  
    [N, Nr, Ns] = shape2_quad(r,s); 
    [J] = jacobian_2d(r,s, x, y, nen); 
    detJ = det(J); 
    invJ = (inv(J));

    B = zeros(nen,2);
    
    Nsr= zeros(2,4); 
    for i = 1:nen
    Nsr(1,i) = Nr(1,i); 
    Nsr(2,i) = Ns(1,i); 
    end
   
    B = (invJ)*(Nsr); 
    B 
   
end

    
    
    

   
    
   
     
        