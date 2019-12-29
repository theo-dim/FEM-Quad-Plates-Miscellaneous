function [Ke] = Ke_heat_tri(D, xn, ien, t, nen, nsd) 
    Ke = zeros(nen); 
    x_bar = 0; 
    y_bar = 0; 
    for i = 1:nen
        x_bar = x_bar + 1/3*xn(1,ien(i)); 
        y_bar = y_bar + 1/3*xn(2,ien(i));
    end;
    
    [N, dNx, dNy, detJ] = shape_tri(x_bar, y_bar, xn, ien, nen); 
    B = B_heat_tri(dNx, dNy, nen, nsd); 
    Ke = (t/2)*detJ*B'*D*B;

    end
