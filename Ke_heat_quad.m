function [Ke] = Ke_heat_quad(D, xn, ien, t, nen, point, weight, nglx, ngly)
    Ke = 0; 
    %need to find r and s
    r = zeros(4,1);
    s = zeros(4,1);
    
    r(1) = point(1,1); 
    r(2) = point(2,1); 
    r(3) = point(2,2);
    r(4) = point(1,2);
    s(1) = point(1,1); 
    s(2) = point(1,2);
    s(3) = point(2,1); 
    s(4) = point(2,2);
    

    
    
    Ketemp = zeros(2); 
    
    for i = 1:4
    [N, B, detJ] = shape_quad(r(i), s(i), xn, ien, nen);
    Ketemp = t*detJ*B'*D*B*weight(i);
    Ketemp
    r(i)
    s(i)
    Ke = Ke+Ketemp; 
    end
    Ke
  
end

  