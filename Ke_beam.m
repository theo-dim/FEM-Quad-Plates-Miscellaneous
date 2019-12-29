%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ke_beam.m - November, 21 2002                           %
% author: Jean H. Prevost                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [Ke,ke,Qe]=Ke_beam(E,A,I,k,xn,ien,nen,ndf,nsd)
%------------------------------------------------------------------------
%  Purpose:
%     compute beam element stiffness
% 
%  Synopsis:
%     [Ke,ke,Qe]=Ke_beam(E,A,xn,ien,nen,ndf,nsd)
% 
%  Variable Description:
%     nen  - number of nodes per element
%     ndf  - number of equations per node
%     nsd  - number of spatial dimensions
%------------------------------------------------------------------------
   node1=ien(1);    node2=ien(2);
     
% compute the length for each element
    v=xn(:,node2)-xn(:,node1)  
    L=norm(v)
 
% local stiffness
    %ke=zeros(2,2);
    ke=zeros(6,6);
    %ke(:,:)=(E*A/L)*[[1 -1]; [-1 1]]; % re-code to include general ke4x4
    format long 
    ke(:,:)=[[E*A/L 0 0 -E*A/L 0 0]; 
            [0 12*E*I/(L^3) 6*E*I/(L^2) 0 -12*E*I/(L^3) 6*E*I/(L^2)];
            [0 6*E*I/(L^2) 4*E*I/L 0 -6*E*I/(L^2) 2*E*I/L]; 
            [-E*A/L 0 0 E*A/L 0 0];
            [0 -12*E*I/(L^3) -6*E*I/(L^2) 0 12*E*I/(L^3) -6*E*I/(L^2)]; 
            [0 6*E*I/(L^2) 2*E*I/L 0 -6*E*I/(L^2) 4*E*I/L]]
% Rotation matrix to global coordinate system
if ( nsd == 1)      % 1D case
    if (xn(1,node1) < xn(1,node2))Qe=1;
    else Qe=-1;
    end;
else

    Qe=zeros(ndf*nen,ndf*nen);
    v=v/norm(v)            % normalize vector
        if (nsd ==2)        % 2D case
            Qe(1,1)=v(1);
            Qe(1,2)=v(2);
            Qe(2,3)=v(1);
            Qe(2,4)=v(2);
        else if (nsd ==3)    % 3D case
            Qe(1,1)=v(1);
            Qe(1,2)=v(2);
            Qe(2,1)=-v(2);
            Qe(2,2)=v(1);
            Qe(3,3)=1;
            Qe(4,4)=v(1);
            Qe(4,5)=v(2);
            Qe(5,4)=-v(2);
            Qe(5,5)=v(1);
            Qe(6,6)=1;
            Qe
        end;
end;
 
% Global Stiffness Ke(i,j)
Ke=zeros(ndf*nen,ndf*nen);
if (nsd > 1)
    Ke(:,:)=Qe(:,:)'*ke(:,:)*Qe(:,:)
else
    Ke(:,:)=Qe*ke(:,:);
end;
 


