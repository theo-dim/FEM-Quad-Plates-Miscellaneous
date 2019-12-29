%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B_heat_quad.m - 1/12/16                                  %
% author: Tehila Stone | Theo Dimitrasopoulos              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B] = B_heat_tri(dNx, dNy, nen, nsd)
    B = zeros(nsd,nen); 
    for i = 1:nen
            B(1,i) = dNx(i); 
    end
    for i = 1:nen
        B(2,i) = dNy(i);
    end
end

        