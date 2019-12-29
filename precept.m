%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% truss.m - October 15 2002                                % 
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
E=1.e8;      % Young's modulus
%%%%%%%%%%%%%
% Geometric %
%%%%%%%%%%%%%
A=0.005;       % Area

%%%%%%%%
% Mesh %
%%%%%%%%
nsd=2;   % number of space dimensions
ndf=nsd; % number of freedom per node 
nen=2;   % number of element nodes

nel=5;   % number of elements/trusses
nnp=4;   % number of nodal points

%%%%%%%%%%%%%%%%%%%%%
% Nodal coordinates %
%%%%%%%%%%%%%%%%%%%%%
% xn(i,N):= coordinate i for node N
% N=1,...,nnp
% i=1,...,nsd
xn=zeros(nsd,nnp); % define coordinates of the system
xn(1,2)=5.0;
xn(2,2)=7.0;
xn(1,3)=11.;
xn(2,3)=7.;
xn(1,4)=16.;

%%%%%%%%%%%%%%%%
% Connectivity %
%%%%%%%%%%%%%%%%
% ien(a,e)=N
% N: global node number - N=1,...,nnp
% e: element number - e=1,...,nel
% a: local node number - a=1,...,nen
ien=zeros(nen,nel);
ien(1,1)=1;  ien(2,1)=2;
ien(1,2)=2;  ien(2,2)=3;
ien(1,3)=3;  ien(2,3)=4;
ien(1,4)=1;  ien(2,4)=3;
ien(1,5)=4;  ien(2,5)=2;

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
idb(1,4)=1;
idb(2,4)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal displacement boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(i,N): prescribed displacement for the dof i of node N
% initialize g
g=zeros(ndf,nnp);
% enter the values
g(1,1)=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal forces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(i,N): prescribed force for the dof i of node N
% initialize f
f=zeros(ndf,nnp);
% enter the values
f(2,2)=-100.;
f(2,3)=-100.;

%---------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number the equations; build the id table %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[id,neq]=number_eq(idb,nnp,ndf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the elemental quantities in the elemental coordinate system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e=1:nel
  [Ke(:,:,e),ke(:,:,e),Qe(:,:,e)]=Ke_truss(E,A,xn,ien(:,e),nen,ndf,nsd);
end;

% Contribution of the prescribed displacements to the elemental force vector
% fe=fe-Ke*Ue
fe=zeros(ndf*nen,nel);     % fe may be non zero in general
Ue=zeros(ndf*nen,nel);
for e=1:nel
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e)=g(i,ien(n,e));
        end
    end
    fe(:,e)=fe(:,e)-Ke(:,:,e)*Ue(:,e);
end

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
        K = addstiff(K,id,Ke(:,:,e),ien(:,e),nen,ndf);
        F = addforce(F,id,fe(:,e),ien(:,e),nen,ndf);
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

% Contribution of the displacement to the elemental force vector
% fe=Ke*Ue
for e=1:nel
    Ue(:,e)=zeros(ndf*nen,1);
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
        end
    end
    fe(:,e)=Ke(:,:,e)*Ue(:,e);
end

% compute reactions R %
R=zeros(ineq,1);
for e=1:nel
    R = addforce(R,idb,fe(:,e),ien(:,e),nen,ndf);
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
    Ue(:,e)=zeros(ndf*nen,1);
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
        end
    end
    if (nsd > 1)
        axial(:,e)=ke(:,:,e)*Qe(:,:,e)*Ue(:,e);
    else
        axial(:,e)=ke(:,:,e)*Qe(e)*Ue(:,e);
    end;
    stress(e)=axial(2,e)/A;
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
plot_results('truss',xn,f,idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp);
