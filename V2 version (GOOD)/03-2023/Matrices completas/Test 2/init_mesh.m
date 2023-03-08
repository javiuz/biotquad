function [x,y]=init_mesh(mesh)

global NN

N=NN;

%N=16;

refinement=2;
x=zeros(N+1,N+1);
y=zeros(N+1,N+1);

if mesh==0 % square regular mesh
    for j=1:N+1
        for i=1:N+1
            x(i,j)=(i-1)*1./N;
            y(i,j)=(j-1)*1./N;
        end
    end

elseif mesh==1 % Random mesh
    x(1,:)=0.;
    x(N+1,:)=1.;
    y(:,1)=0.;
    y(:,N+1)=1.;

    % Artículo COMG_2020
    
    for j=2:N
        r=rand;
        y(1,j)=(j-1.3)*1./N+0.6*r*1./N;
        r=rand;
        y(N+1,j)=(j-1.3)*1./N+0.6*r*1./N;
    end

    for i=2:N
        r=rand;
        x(i,1)=(i-1.3)*1./N+0.6*r*1./N;
        r=rand;
        x(i,N+1)=(i-1.3)*1./N+0.6*r*1./N;
    end

    for j=2:N
        for i=2:N
            r=rand;
            x(i,j)=(i-1.3)*1./N+0.6*r*1./N;
            r=rand;
            y(i,j)=(j-1.3)*1./N+0.6*r*1./N;
        end 
    end
    
    % Usado en mi TFM

%     for j=2:N
%         r=rand;
%         y(1,j)=(j-5./4.)*1./N+sqrt(2.)/(3.*N)*r;
%         r=rand;
%         y(N+1,j)=(j-5./4.)*1./N+sqrt(2.)/(3.*N)*r;
%     end
% 
%     for i=2:N
%         r=rand;
%         x(i,1)=(i-5./4.)*1./N+sqrt(2.)/(3.*N)*r;
%         r=rand;
%         x(i,N+1)=(i-5./4.)*1./N+sqrt(2.)/(3.*N)*r;
%     end
% 
%     for j=2:N
%         for i=2:N
%             r=rand;
%             x(i,j)=(i-5./4.)*1./N+sqrt(2.)/(3.*N)*r;
%             r=rand;
%             y(i,j)=(j-5./4.)*1./N+sqrt(2.)/(3.*N)*r;
%         end 
%     end
    
elseif mesh==11 % Random mesh 3
    x(1,:)=0.;
    x(N+1,:)=1.;
    y(:,1)=0.;
    y(:,N+1)=1.;

%    call random_seed()

    for j=2:N
        r=rand;
        y(1,j)=(j-4./3.)*1./N+2./(3.*N)*r;
        r=rand;
        y(N+1,j)=(j-4./3.)*1./N+2./(3.*N)*r;
    end

    for i=2:N
        r=rand;
        x(i,1)=(i-4./3.)*1./N+2./(3.*N)*r;
        r=rand;
        x(i,N+1)=(i-4./3.)*1./N+2./(3.*N)*r;
    end

    for j=2:N
        for i=2:N
            r=rand;
            x(i,j)=(i-4./3.)*1./N+2./(3.*N)*r;
            r=rand;
            y(i,j)=(j-4./3.)*1./N+2./(3.*N)*r;
        end 
    end
elseif mesh==2 % D. Boffi mesh
    hxx=1./N;
    hyy=1./N;
    hysmall=1./(2.*N);
    hybig=3./(2.*N);
    
    for i=1:N+1
        x(i,1:N+1)=(i-1)*hxx;
    end 
    for j=1:2:N+1
        y(1:N+1,j)=(j-1)*hyy;
    end
    for j=2:2:N
        for i=1:2:N+1
            y(i,j)=y(i,j-1)+hysmall;
        end
        for i=2:2:N
            y(i,j)=y(i,j-1)+hybig;
        end
    end

elseif mesh==3 % Regular mesh O(h^2) parallelograms
    x(1,:)=0.;
    x(N+1,:)=1.;
    y(:,1)=0.;
    y(:,N+1)=1.;
    
    c1=0.06;
    c2=-0.05;
    
    for j=2:N
        y(1,j)=(j-1)*1./N;
        y(N+1,j)=(j-1)*1./N;
    end
    
    for i=2:N
        x(i,1)=(i-1)*1./N;
        x(i,N+1)=(i-1)*1./N;
    end

    for j=2:N
        for i=2:N
            xx=(i-1)*1./N;
            yy=(j-1)*1./N;
            x(i,j)=xx+c1*sin(2.*pi*xx)*sin(2.*pi*yy);
            y(i,j)=yy+c2*sin(2.*pi*xx)*sin(2.*pi*yy);
        end
    end
% else % Porous media example mesh
%     for j=1:N+1
%         for i=1:N+1
%             x(i,j)=(i-1)*1./N;
%         end 
%     end 
%     for i=1:N+1
%         ys=sqrt(1.21-((i-1)*1./N-0.1)^2)-0.4;
%         yn=sqrt(1.44-((i-1)*1./N-0.1)^2)-0.4;
%         hs=ys/(9.*mm);
%         hc=(yn-ys)/(2.*mm);
%         hn=(1.-yn)/(5.*mm);
%         for j=1:9*mm+1
%             y(i,j)=(j-1)*hs;
%         end
%         for j=9*mm+2:11*mm+1
%             y(i,j)=ys+(j-9*mm-1)*hc;
%         end
%         for j=11*mm+2:16*mm+1
%             y(i,j)=yn+(j-11*mm-1)*hn;
%         end
%     end
% end

elseif mesh==4  % Kershaw mesh
                    % N must be multiple of 8.
 for j=1:N+1
    for i=1:N+1
        y(i,j)=(j-1)/N;
    end
end

d=0.4;

for j=1:N/8+1
    hxl=(1.-d)*2/N;
    hxr=d*2./N;
    for i=1:N/2+1
        x(i,j)=(i-1)*hxl;
    end
    for i=N/2+2:N+1
        x(i,j)=(1.-d)+(i-N/2-1)*hxr;
    end
end
j=N/4+1;
hxl=d*2./N;
hxr=(1.-d)*2/N;
for i=1:N/2+1
    x(i,j)=(i-1)*hxl;
end
for i=N/2+2:N+1
    x(i,j)=d+(i-N/2-1)*hxr;
end
j=3*N/4+1;
hxl=(1.-d)*2./N;
hxr=d*2/N;
for i=1:N/2+1
    x(i,j)=(i-1)*hxl;
end
for i=N/2+2:N+1
    x(i,j)=1.-d+(i-N/2-1)*hxr;
end
for j=7*N/8+1:N+1
    hxl=d*2./N;
    hxr=(1.-d)*2/N;
    for i=1:N/2+1
        x(i,j)=(i-1)*hxl;
    end
    for i=N/2+2:N+1
        x(i,j)=d+(i-N/2-1)*hxr;
    end
end

for i=1:N+1
    ds=x(i,N/4+1)-x(i,N/8+1);
    dc=x(i,3*N/4+1)-x(i,N/4+1);
    dn=x(i,7*N/8+1)-x(i,3*N/4+1);
    for j=N/8+2:N/4
        x(i,j)=x(i,N/8+1)+(j-N/8-1)*ds/(N/8);
    end
    for j=N/4+2:3*N/4
        x(i,j)=x(i,N/4+1)+(j-N/4-1)*dc/(N/2);
    end
    for j=3*N/4+2:7*N/8
        x(i,j)=x(i,3*N/4+1)+(j-3*N/4-1)*dn/(N/8);
    end
end

elseif mesh==21 %FVCA5_mesh1
    
    A = load('mesh_1.m');
    [n1 n2] = size(A);
    n = sqrt(n1);
    ind = 1;
    for j = 1:n
        for i = 1:n
            x(i,j) = A(ind,1);
            y(i,j) = A(ind,2);
            ind = ind+1;
        end
    end
%     figure
%     hold on
%     for j=1:n
%         plot(x(:,j),y(:,j),'k')
%     end
%     for i=1:n
%         plot(x(i,:),y(i,:),'k')
%     end
%     hold off
 
    nf = n;
    for k = 1:refinement
        nf = 2*nf-1;
        for j = 1:2:nf
            for i = 1:2:nf
                xf(i,j) = x((i+1)/2,(j+1)/2);
                yf(i,j) = y((i+1)/2,(j+1)/2);
            end
            for i = 2:2:nf
                xf(i,j) = 0.5*(xf(i-1,j)+xf(i+1,j));
                yf(i,j) = 0.5*(yf(i-1,j)+yf(i+1,j));             
            end
        end
        for j = 2:2:nf
            for i = 1:nf
                xf(i,j) = 0.5*(xf(i,j-1)+xf(i,j+1));
                yf(i,j) = 0.5*(yf(i,j-1)+yf(i,j+1));
            end
        end  
        x = xf;
        y = yf;
    end
%     figure
%     hold on
%     for j=1:nf
%         plot(x(:,j),y(:,j),'k')
%     end
%     for i=1:nf
%         plot(x(i,:),y(i,:),'k')
%     end
%     hold off   
    
elseif mesh==22 %FVCA5_mesh2
    A = load('mesh_2.m');
    [n1 n2] = size(A);
    n = sqrt(n1);
    ind = 1;
    for j = 1:n
        for i = 1:n
            x(i,j) = A(ind,1);
            y(i,j) = A(ind,2);
            ind = ind+1;
        end
    end
%     figure
%     hold on
%     for j=1:n
%         plot(x(:,j),y(:,j),'k')
%     end
%     for i=1:n
%         plot(x(i,:),y(i,:),'k')
%     end
%     hold off
 
    nf = n;
    for k = 1:refinement
        nf = 2*nf-1;
        for j = 1:2:nf
            for i = 1:2:nf
                xf(i,j) = x((i+1)/2,(j+1)/2);
                yf(i,j) = y((i+1)/2,(j+1)/2);
            end
            for i = 2:2:nf
                xf(i,j) = 0.5*(xf(i-1,j)+xf(i+1,j));
                yf(i,j) = 0.5*(yf(i-1,j)+yf(i+1,j));             
            end
        end
        for j = 2:2:nf
            for i = 1:nf
                xf(i,j) = 0.5*(xf(i,j-1)+xf(i,j+1));
                yf(i,j) = 0.5*(yf(i,j-1)+yf(i,j+1));
            end
        end  
        x = xf;
        y = yf;
    end
    
end

% Dibujo de la malla

% fig=figure;
%     hold on
%     for j=1:N+1
%         plot(x(:,j),y(:,j),'k','LineWidth',1.2)
%     end
%     for i=1:N+1
%         plot(x(i,:),y(i,:),'k','LineWidth',1.2)
%     end
%     set(gca, 'YTick', []);
% set(gca, 'XTick', []);
% set(gca,'PlotBoxAspectRatio',[1 1 1]);
% %pos = get(gcf,'paperposition');
% 
% %set(gcf,'paperposition',[pos(1),pos(2), 4, 4]);
% orient(fig,'landscape');
% print(fig,'prueba.pdf','-dpdf')
% 
%     hold off 

return
end