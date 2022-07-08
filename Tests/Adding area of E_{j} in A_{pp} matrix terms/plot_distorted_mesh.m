N=4;

xx=zeros(N+1,N+1);
yy=zeros(N+1,N+1);

for j=1:N+1
        for i=1:N+1
            xx(i,j)=(i-1)*1./N;
            yy(i,j)=(j-1)*1./N;
        end
end

x=xx+0.03*cos(3*pi*xx).*cos(3*pi*yy);
y=yy-0.04*cos(3*pi*xx).*cos(3*pi*yy);

hold on
for j=1:N+1
      plot(x(:,j),y(:,j),'k','LineWidth',1.2)
end
for i=1:N+1
      plot(x(i,:),y(i,:),'k','LineWidth',1.2)
end
axis equal

refinement=2;

nf = N;
    for k = 1:refinement
        nf = 2*nf;
        for j = 1:2:nf+1
            for i = 1:2:nf+1
                xf(i,j) = x((i+1)/2,(j+1)/2);
                yf(i,j) = y((i+1)/2,(j+1)/2);
            end
            for i = 2:2:nf+1
                xf(i,j) = 0.5*(xf(i-1,j)+xf(i+1,j));
                yf(i,j) = 0.5*(yf(i-1,j)+yf(i+1,j));             
            end
        end
        for j = 2:2:nf+1
            for i = 1:nf+1
                xf(i,j) = 0.5*(xf(i,j-1)+xf(i,j+1));
                yf(i,j) = 0.5*(yf(i,j-1)+yf(i,j+1));
            end
        end  
        x = xf;
        y = yf;
    end
    
figure
hold on
for j=1:nf+1
      plot(x(:,j),y(:,j),'k','LineWidth',1.2)
end
for i=1:nf+1
      plot(x(i,:),y(i,:),'k','LineWidth',1.2)
end
axis equal
