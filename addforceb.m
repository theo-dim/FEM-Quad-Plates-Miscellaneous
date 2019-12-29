%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addforce.m - November, 21 2002                           %
% author: Jean H. Prevost                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fb]=addforceb(Fb,id,feb,ien,nen,ndf)

%------------------------------------------------------------------------
%  Purpose:
%     add elemental force in global force vector
%
%  Synopsis:
%     [F]=addforce(id,fe,ien,nen,ndf)
%
%  Variable Description:
%     nen  - number of nodes per element
%     ndf  - number of equations per node
%------------------------------------------------------------------------

for n=1:nen
    for i=1:ndf
        if (id(i,ien(n)) > 0)
            P=id(i,ien(n));
            Fb(P)=Fb(P)+feb(i+(n-1)*ndf);
        end
    end
end

