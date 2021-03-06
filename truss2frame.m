%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% truss.m - October 15 2002                                %
%           revised: September 21, 2006                    %
% author: Jean H. Prevost + David Luet                     %
% analyses of 1,2 and 3D elastic trusses                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; % removes all variables from the workspace.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          DATA                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
% Material %
%%%%%%%%%%%%
E=30;      % Young's modulus in kN/m^2

%%%%%%%%%%%%%
% Geometric %
%%%%%%%%%%%%%
A=1500;       % Area in m^2
I = 2*10^9; %Moment of inertia given in mm^4

%%%%%%%%
% Mesh %
%%%%%%%%
nsd=3;   % number of space dimensions
ndf=nsd; % number of freedom per node 
nen=3;   % number of element nodes

nel=3;   % number of elements/trusses
nnp=4;   % number of nodal points

%%%%%%%%%%%%%%%%%%%%%
% Nodal coordinates %
%%%%%%%%%%%%%%%%%%%%%
% xn(i,N):= coordinate i for node N in m
% N=1,...,nnp
% i=1,...,nsd
xn=zeros(nsd,nnp);
xn(1,1)=0;
xn(2,1)=0;
xn(1,2)=0;
xn(2,2)=6.0;
xn(1,3)=12.0;
xn(2,3)=6.0;
xn(1,4)=12;
xn(2,4)=0;

%%%%%%%%%%%%%%%%
% Connectivity %
%%%%%%%%%%%%%%%%
% ien(a,e)=N
% N: global node number - N=1,...,nnp
% e: element number - e=1,...,nel
% a: local node number - a=1,...,nen
% mat(e) is used to identify the area of each truss element
ien=zeros(nen,nel); mat(nel)=0;
ien(1,1)=1;  ien(2,1)=2;  mat(1)=1;
ien(1,2)=2;  ien(2,2)=3;  mat(2)=1;
ien(1,3)=3;  ien(2,3)=4;  mat(2)=1;

%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%
% prescribed displacement (essential boundary condition)
%
% idb(i,N)=1 if the degree of freedom i of the node N is prescribed
%         =0 otherwise
%
% 1) initialize idb to 0
idb=zeros(ndf,nnp);
% 2) enter the flag for prescribed displacement boundary conditions
idb(1,1)=1;
idb(2,1)=1;
idb(1,2)=0;
idb(2,2)=0;
idb(1,3)=0;
idb(2,3)=0;
idb(1,4)=1;
idb(2,4)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal displacement boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(i,N): prescribed displacement for the dof i of node N
% initialize g
g=zeros(ndf,nnp);
% enter the values in m

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal forces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(i,N): prescribed force for the dof i of node N
% initialize f
f=zeros(ndf,nnp);
% enter the values in kN
f(2,2)=10;

%---------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number the equations; build the id table %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[id,neq]=number_eq(idb,nnp,ndf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the elemental quantities in the elemental coordinate system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e=1:nel
  [Ke(:,:,e),ke(:,:,e),Qe(:,:,e)]=Ke_truss(E,A(mat(e)),xn,ien(:,e),nen,ndf,nsd);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Assembly operation                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------
% build K and F
%----------------
K=zeros(neq,neq);
F=zeros(neq,1);

% input the prescribed nodal forces in F
for N=1:nnp
    for i=1:ndf
        if (id(i,N) > 0)
            P=id(i,N);
            F(P)=f(i,N);
%            F(id(i,N))=f(i,N);
        end
    end
end

% compute global K and F
if (neq > 0)
    for e=1:nel
      % stiffness contribution
        K = addstiff(K,id,Ke(:,:,e),ien(:,e),nen,ndf);
      % Contribution of the prescribed displacements: fe=fe-Ke*Ue
        fe=zeros(ndf*nen,1);    
        Ue=zeros(ndf*nen,1);
        for n=1:nen
            for i=1:ndf
                Ue(i+(n-1)*ndf)=g(i,ien(n,e));
            end
        end
        fe(:) = fe(:) - Ke(:,:,e)*Ue(:);
      % add to force array
        F = addforce(F,id,fe(:),ien(:,e),nen,ndf);
    end
end

% Solve the system
if (neq > 0)
    U=K\F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     post processing                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% complete U %
%%%%%%%%%%%%%%
Ucomp=zeros(ndf,nnp);
for N=1:nnp
    for i=1:ndf
        if (id(i,N) == 0) 
            Ucomp(i,N)=g(i,N);
        else
            P=id(i,N);
            Ucomp(i,N)=U(P);
        end
    end
end
% print results
disp('Nodal Displacements:')
disp(' node     d1     d2')
for N=1:nnp
   disp(sprintf('%5d %7g %7g',N,Ucomp(:,N)))
end 
disp(' ')

%%%%%%%%%%%%%
% REACTIONS %
%%%%%%%%%%%%%
% build the idb table; overwrite original idb table
% idb(i,N): equation number associated with dof i of node N
ineq=0; % number of equations
for i=1:ndf
    for N=1:nnp
        if (idb(i,N) > 0)  % assign an equation number to all prescribed nodes
            ineq=ineq+1;
            idb(i,N)=ineq;
        end;
    end;
end;

% compute reactions R %
R=zeros(ineq,1);
for e=1:nel
  % Contribution of the displacement:  fe=Ke*Ue
    Ue=zeros(ndf*nen,1);
    fe=zeros(ndf*nen,1); 
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf)=Ucomp(i,ien(n,e));
        end
    end
    fe(:) = fe(:) + Ke(:,:,e)*Ue(:);
    R = addforce(R,idb,fe(:),ien(:,e),nen,ndf);
end

% Collect reactions
Rcomp=zeros(ndf,nnp);
for N=1:nnp
    for i=1:ndf
        if (idb(i,N) > 0)
            Rcomp(i,N)=R(idb(i,N));
        end
    end
end
% print results
disp('Nodal Reactions')
disp(' node     R1     R2')
for N=1:nnp
   disp(sprintf('%5d %7g %7g',N,Rcomp(:,N)))
end 
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%
% AXIAL FORCES/STRESSES %
%%%%%%%%%%%%%%%%%%%%%%%%%
for e=1:nel
    Ue=zeros(ndf*nen,1);
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf)=Ucomp(i,ien(n,e));
        end
    end
    if (nsd > 1)
        axial(:,e)=ke(:,:,e)*Qe(:,:,e)*Ue(:);
    else
        axial(:,e)=ke(:,:,e)*Qe(e)*Ue(:);
    end;
    stress(e)=axial(2,e)/A(mat(e));
    strain(e)=stress(e)/E;
end;
% print results
disp('Element Axial force/stress/strain')
disp(' elem  force  stress  strain')
for e=1:nel
   disp(sprintf('%5d %7g %7g %7g',e,axial(2,e),stress(e),strain(e)))
end 
disp(' ')

%%%%%%%%%%%%%%%%%%%%
% plot the results %
%%%%%%%%%%%%%%%%%%%%
plot_results('beam',xn,f,idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp);
