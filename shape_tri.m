function [N,dNx, dNy, detJ] = shape_tri(x_bar, y_bar, xn, ien, nen)

N= zeros(1, nen);
dNx = zeros(1,nen);
dNy = zeros(1,nen);
x = zeros(1,nen);
y = zeros(1,nen);
    for i= 1:nen
        x(i) = xn(1, ien(i)); 
        y(i) = xn(2, ien(i)); 
    end
    A =  (1/2)*(x(2)*y(3)+x(1)*y(2)+x(3)*y(1)-x(2)*y(1)-x(1)*y(3)-x(3)*y(2)); 
    % N1 = (1/(2*A))*((x(2)*y(3)-x(3)*y(2))+(y(2)-y(3))*x_bar+(x(3)-x(2))*y_bar);
    % N2 = (1/(2*A))*((x(3)*y(1)-x(1)*y(3))+(y(3)-y(1))*x_bar+(x(1)-x(3))*y_bar);
    % N3 = (1/(2*A))*((x(1)*y(2)-x(2)*y(1))+(y(1)-y(2))*x_bar+(x(2)-x(1))*y_bar);
    
    N1x = (1/(2*A))*(y(2)-y(3)); %y(2)-y(1)
    N2x = (1/(2*A))*(y(3)-y(1));
    N3x = (1/(2*A))*(y(1)-y(2));
    N1y = (1/(2*A))*(x(3)-x(2));
    N2y = (1/(2*A))*(x(1)-x(3));
    N3y = (1/(2*A))*(x(2)-x(1));
    
    detJ =2*A;
    
    dNx(1) = N1x;
    dNx(2) = N2x;
    dNx(3) = N3x;
    dNy(1) = N1y;
    dNy(2) = N2y;
    dNy(3) = N3y;
        
end

