%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_results_elas.m - November,23 2002                    % 
% author: David Luet - luet@princeton.edu                   %
% plots results and b.c. for trusses,beams, heat conduction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_results_elas(type,xn,f,Idb,Ucomp,Rcomp,sig,epsilon,sig1,ien,nel,nen,nsd,ndf,nnp);

display_factor=0.1;
eps=1e-4;

% characteristic distances
xmax=max(max(xn(1,:)));
xmin=min(min(xn(1,:)));

ymax=max(max(xn(2,:)));
ymin=min(min(xn(2,:)));

% Lcar=max([xmax-xmin;ymax-ymin]);
Lcar=((xmax-xmin)+(ymax-ymin))/2;

k=0;
while (k~=14)
    k=menu('results','undeformed mesh and BCs', 'undeformed and deformed mesh',...
        'reactions','sigxx','sigyy', 'tauxy', ...
        'epsxx','epsyy', 'gamxy', ...
        'sig1','sig2','tau','theta','exit');
    
    switch k
    case 1
        figure;    
        axis equal;
        title('Undeformed mesh and BCs');
        hold on;
        plot_mesh_underformed(nel,ien,xn,nnp,nsd,nen);
        numbers(nel,ien,xn,nnp,nsd,nen);
        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd);
        plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd);
        hold off;
    case 2
        figure;
        axis equal;
        title('Undeformed and deformed mesh');
        hold on;
        plot_mesh_underformed(nel,ien,xn,nnp,nsd,nen);
        plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,nnp);
        hold off;
    case 3
        figure;
        axis equal;
        title('reactions');
        hold on;
        plot_mesh_underformed(nel,ien,xn,nnp,nsd,nen);
        plot_reactions(type,Lcar,Rcomp,Idb,display_factor,xn,nnp,ndf,nsd)
        hold off;
    case 4
        figure;
        axis equal;
        title('\sigma_{xx}');
        hold on;
        plot_stress_strain(sig, 1, nel,nen, nsd, ien, xn);
        hold off;
    case 5
        figure;
        axis equal;
        title('\sigma_{yy}');
        hold on;
        plot_stress_strain(sig, 2, nel,nen, nsd, ien, xn);
        hold off;
    case 6
        figure;
        axis equal;
        title('\tau_{xy}');
        hold on;
        plot_stress_strain(sig, 3, nel,nen, nsd, ien, xn);
        hold off;
    case 7
        figure;
        axis equal;
        title('\epsilon_{xx}');
        hold on;
        plot_stress_strain(epsilon, 1, nel,nen, nsd, ien, xn);
        hold off;
    case 8
        figure;
        axis equal;
        title('\epsilon_{yy}');
        hold on;
        plot_stress_strain(epsilon, 2, nel,nen, nsd, ien, xn);
        hold off;
    case 9
        figure;
        axis equal;
        title('\gamma_{xy}');
        hold on;
        plot_stress_strain(epsilon, 3, nel,nen, nsd, ien, xn);
        hold off;
    case 10
        figure;
        axis equal;
        title('\sigma_I');
        hold on;
        plot_stress_strain(sig1, 1, nel,nen, nsd, ien, xn);
        hold off;
    case 11
        figure;
        axis equal;
        title('\sigma_{II}');
        hold on;
        plot_stress_strain(sig1, 2, nel,nen, nsd, ien, xn);
        hold off;
    case 12
        figure;
        axis equal;
        title('\tau');
        hold on;
        plot_stress_strain(sig1, 3, nel,nen, nsd, ien, xn);
        hold off;
    case 13
        figure;
        axis equal;
        title('\theta');
        hold on;
        plot_stress_strain(sig1, 4, nel,nen, nsd, ien, xn);
        hold off;
    end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       PLOT FUNCTIONS                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% undeformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mesh_underformed(nel,ien,xn,nnp,nsd,nen)
if (nen == 3)
    triplot(ien',xn(1,:),xn(2,:),'blue');
else
    for e=1:nel,
        for m=1:nen,
            node(m)=ien(m,e);
        end;
        x=[xn(1,node(1)); xn(1,node(2)); xn(1,node(3)); xn(1,node(4)); xn(1,node(1))];
        y=[xn(2,node(1)); xn(2,node(2)); xn(2,node(3)); xn(2,node(4)); xn(2,node(1))];
        plot(x,y,'b-o');
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbers the nodes and elements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numbers(nel,ien,xn,nnp,nsd,nen)
% number the nodes 
for n=1:nnp
    if (nsd > 1)
        text(xn(1,n),xn(2,n), num2str(n),'FontSize',10);
    else
        text(xn(1,n),0, num2str(n),'FontSize',10);
    end;
end
% number the elements 
% assume nsd=2
for e=1:nel        
    xg=zeros(2,1);
    for i=1:nen
        xg=xg+xn(:,ien(i,e));
    end
    xg=xg/nen;
    s=sprintf('(%d)', e);
    text(xg(1),xg(2), s,'FontSize',10);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% deformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,nnp);
scale=display_factor*Lcar/max(max(abs(Ucomp))); % scale factor for the displacements
switch nen
case 3
    for i=1:nnp,
        xt(:,i)=xn(:,i)+scale*Ucomp(:,i);
    end
    triplot(ien',xt(1,:),xt(2,:),'red');
case 4
    for e=1:nel
        for i=1:4
            node(i)=ien(i,e);
            xt(:,i)=xn(:,node(i))+scale*Ucomp(:,node(i));
        end;
        plot(xt(1,:),xt(2,:)','r-o');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions on displacements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd)
alpha=display_factor*Lcar; % scale factor  for bc symbols
for P=1:nnp
    if (nsd >1)
        if ((Idb(1,P) ~= 0) & (Idb(2,P) ~= 0))
            bc_symbols(xn(:,P),alpha,3);
        end;
        if ((Idb(1,P) ~= 0) & (Idb(2,P) == 0))
            bc_symbols(xn(:,P),alpha,2);
        end;
        if ((Idb(1,P) == 0) & (Idb(2,P) ~= 0))
            bc_symbols(xn(:,P),alpha,1);
        end;
    else
        if (Idb(1,P) ~= 0) bc_symbols([xn(1,P),0],alpha,2); 
        end;
    end;       
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions on force %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd)
delta=display_factor*Lcar/max(max(abs(f))); % scale factor for force b.c.
alpha=2*display_factor*Lcar; % scale factor for moment
for N=1:nnp
    if ( nsd >1)
        if ( (f(1,N) ~=0 ) | (f(2,N) ~= 0))
            quiver(xn(1,N),xn(2,N),f(1,N),f(2,N),delta,'r');
        end;        
    else
        if (f(1,N) ~=0) quiver(xn(1,N),0,f(1,N),0,delta,'r');
        end;
    end;
end;



%%%%%%%%%%%%%
% reactions %
%%%%%%%%%%%%%
function plot_reactions(type,Lcar,Rcomp,Idb,display_factor,xn,nnp,ndf,nsd)
beta=display_factor*Lcar/max(max(abs(Rcomp))); % scale factor for reactions
alpha=2*display_factor*Lcar; % scale factor for moment
for N=1:nnp
    RN=zeros(ndf);
    if ((Idb(1,N) ~= 0) | (Idb(2,N) ~= 0))
        quiver(xn(1,N),xn(2,N),Rcomp(1,N),Rcomp(2,N),beta,'r');
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%
% strains and stresses %
%%%%%%%%%%%%%%%%%%%%%%%%
function plot_stress_strain(tab, ind, nel, nen, nsd, ien, xn)
m=nen+1;
nsize=3; % size of the stress and strain vectors
maxi=max(tab(ind,:));
mini=min(tab(ind,:));
mean=0.5*(maxi+mini);
span=maxi-mini;
if ((abs(span) < 0.05*abs(mean)))
    uniform=1;
else
    uniform=0;
end;

if (nen == 3)
    for e=1:nel
        for i=1:nen
            node(i)=ien(i,e);
        end
        x=[xn(1,node(1)); xn(1,node(2));xn(1,node(3));xn(1,node(1))];
        y=[xn(2,node(1)); xn(2,node(2));xn(2,node(3));xn(2,node(1))];
        if (uniform ==0)
            c=[tab(ind,e); tab(ind,e); tab(ind,e);tab(ind,e)];
        else
            c=[mean; mean; mean; mean];
        end;
        fill(x,y,c);
    end
    colorbar; 
elseif (nen == 4)
    for e=1:nel
        for i=1:nen
            node(i)=ien(i,e);
        end
        x=[xn(1,node(1)); xn(1,node(2));xn(1,node(3));xn(1,node(4));xn(1,node(1))];
        y=[xn(2,node(1)); xn(2,node(2));xn(2,node(3));xn(2,node(4));xn(2,node(1))];
        if (uniform == 0)
            c=[tab(ind,e); tab(ind,e); tab(ind,e);tab(ind,e); tab(ind,e)];
        else
            c=[mean;mean;mean;mean;mean];
        end
        fill(x,y,c);
    end
    colorbar; 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              BOUNDARY CONDITIONS SYMBOLS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bc_symbols(xp,alpha,symbol);

switch symbol
case 1
    % v fixed
    x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
    y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
    
    line(x,y,'Color','k','LineWidth',1.2);
    
    for i=0:3,
        circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(7/8)*alpha],alpha/8);
    end;
    
    
case 2
    % u fixed
    x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
    y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];
    
    line(x,y,'Color','k','LineWidth',1.2);
    
    for i=0:3,
        circle([xp(1)-(7/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
    end;
    
case 3  
    % u and v fixed
    x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
    y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
    
    line(x,y,'Color','k','LineWidth',1.2);
    
    for i=0:3,
        line([xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)], ...
            [xp(2)-(3/4)*alpha;xp(2)-alpha],'Color','k','LineWidth',1.2);
    end;
    
case 4  
    % v and theta fixed 
    x=[xp(1)-alpha/2;xp(1)+alpha/2];
    y=[xp(2);xp(2)];
    
    line(x,y,'Color','k','LineWidth',1.2);
    
    for i=0:3,
        circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(1/8)*alpha],alpha/8);
    end;
    
case 5
    % u and theta fixed
    x=[xp(1);xp(1)];
    y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];
    
    line(x,y,'Color','k','LineWidth',1.2);
    
    for i=0:3,
        circle([xp(1)-(1/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
    end;    
    
case 6
    % u, v and theta fixed
    line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],'Color','k','LineWidth',1.2);
    for i=0:3,
        line([xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4],[xp(2),xp(2)-alpha/4]...
            ,'Color','k','LineWidth',1.2);
    end;
    
end;
   
    

function circle(x0,r);
theta=0:0.1:2*pi;
x=r*cos(theta)+x0(1);
y=r*sin(theta)+x0(2);

plot(x,y,'k','LineWidth',1.2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       MOMENTS                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_moments(xp,alpha,sign);
switch sign
case 0  %positive moment
    r=alpha/4;
    theta=0:0.1:3*pi/2;
    x=r*cos(theta)+xp(1);
    y=r*sin(theta)+xp(2);
    
    plot(x,y,'k','LineWidth',1.2);
    line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
        'Color','k','LineWidth',1.2);
    line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
        'Color','k','LineWidth',1.2);

case 1  % negative moment
    r=alpha/4;
    theta=pi:-0.1:-pi/2;
    x=r*cos(theta)+xp(1);
    y=r*sin(theta)+xp(2);
    
    plot(x,y,'k','LineWidth',1.2);
    line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
        'Color','k','LineWidth',1.2);
    line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
        'Color','k','LineWidth',1.2);
end;