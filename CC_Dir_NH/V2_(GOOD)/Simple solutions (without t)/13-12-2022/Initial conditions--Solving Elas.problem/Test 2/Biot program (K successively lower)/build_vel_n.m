function vel_n=build_vel_n(z)

global NN 

vel_n=zeros(NN,NN,8);

ind=1; % denota el índice para movernos en las componentes del vector z

% SW corner node

i=1;
j=1;
vel_n(i,j,1)=z(ind,1);
vel_n(i,j,2)=z(ind+1,1);
ind=ind+2;

% South nodes (j=1)

for i=2:NN
    vel_n(i,j,1)=z(ind,1);
    vel_n(i,j,2)=z(ind+1,1);
    vel_n(i-1,j,3)=z(ind+1,1);
    vel_n(i-1,j,4)=z(ind+2,1);
    ind=ind+3;
end

% SE corner node (j=1)

i=NN+1;
vel_n(i-1,j,3)=z(ind,1);
vel_n(i-1,j,4)=z(ind+1,1);
ind=ind+2;

for j=2:NN
    
    % West nodes 
    
    i=1;
    vel_n(i,j-1,8)=z(ind,1);
    vel_n(i,j-1,7)=z(ind+1,1);
    vel_n(i,j,1)=z(ind+1,1);
    vel_n(i,j,2)=z(ind+2,1);
    ind=ind+3;

    % Central nodes
    
    for i=2:NN
        vel_n(i-1,j-1,5)=z(ind,1);
        vel_n(i,j-1,8)=z(ind,1);
        vel_n(i,j-1,7)=z(ind+1,1);
        vel_n(i,j,1)=z(ind+1,1);
        vel_n(i,j,2)=z(ind+2,1);
        vel_n(i-1,j,3)=z(ind+2,1);
        vel_n(i-1,j,4)=z(ind+3,1);
        vel_n(i-1,j-1,6)=z(ind+3,1);
        ind=ind+4;
    end
    
    % East nodes 
    
    i=NN+1;
    vel_n(i-1,j-1,5)=z(ind,1);   
    vel_n(i-1,j,3)=z(ind+1,1);  
    vel_n(i-1,j,4)=z(ind+2,1); 
    vel_n(i-1,j-1,6)=z(ind+2,1);
    ind=ind+3;
end

% NW corner node

i=1;
j=NN+1;

vel_n(i,j-1,8)=z(ind,1);
vel_n(i,j-1,7)=z(ind+1,1);
ind=ind+2;

% North nodes (j=N+1)

for i=2:NN
    vel_n(i-1,j-1,5)=z(ind,1);   
    vel_n(i,j-1,8)=z(ind,1);   
    vel_n(i,j-1,7)=z(ind+1,1);  
    vel_n(i-1,j-1,6)=z(ind+2,1);
    ind=ind+3;
end

% NE corner node

i=NN+1;
vel_n(i-1,j-1,5)=z(ind,1);
vel_n(i-1,j-1,6)=z(ind+1,1);

return
end

% function l=side_length(i,j,k,l)
% global x y
% l=sqrt((x(i,j)-x(k,l))^2+(y(i,j)-y(k,l))^2);
% return
% end