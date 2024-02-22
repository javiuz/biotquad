function sigma_n=build_sigma_n_2(sigma)

global NN 

sigma_n=zeros(NN,NN,16);

ind=1; % denota el índice para movernos en las componentes del vector sigma

% SW corner node

i=1;
j=1;
sigma_n(i,j,1)=sigma(ind,1);
sigma_n(i,j,2)=sigma(ind+1,1);
sigma_n(i,j,3)=sigma(ind+2,1);
sigma_n(i,j,4)=sigma(ind+3,1);
ind=ind+4;

% South nodes (j=1)

for i=2:NN
    sigma_n(i,j,1)=sigma(ind,1);
    sigma_n(i,j,2)=sigma(ind+1,1);
    sigma_n(i,j,3)=sigma(ind+2,1);
    sigma_n(i,j,4)=sigma(ind+3,1);
    sigma_n(i-1,j,5)=sigma(ind+2,1);
    sigma_n(i-1,j,6)=sigma(ind+3,1);
    sigma_n(i-1,j,7)=sigma(ind+4,1);
    sigma_n(i-1,j,8)=sigma(ind+5,1);
    ind=ind+6;
end

% SE corner node (j=1)

i=NN+1;
sigma_n(i-1,j,5)=sigma(ind,1);
sigma_n(i-1,j,6)=sigma(ind+1,1);
sigma_n(i-1,j,7)=sigma(ind+2,1);
sigma_n(i-1,j,8)=sigma(ind+3,1);
ind=ind+4;

for j=2:NN
    
    % West nodes 
    
    i=1;
    sigma_n(i,j-1,15)=sigma(ind,1);
    sigma_n(i,j-1,16)=sigma(ind+1,1);
    sigma_n(i,j-1,13)=sigma(ind+2,1);
    sigma_n(i,j-1,14)=sigma(ind+3,1);
    sigma_n(i,j,1)=sigma(ind+2,1);
    sigma_n(i,j,2)=sigma(ind+3,1);
    sigma_n(i,j,3)=sigma(ind+4,1);
    sigma_n(i,j,4)=sigma(ind+5,1);
    ind=ind+6;

    % Central nodes
    
    for i=2:NN
        sigma_n(i-1,j-1,9)=sigma(ind,1);
        sigma_n(i-1,j-1,10)=sigma(ind+1,1);
        sigma_n(i,j-1,15)=sigma(ind,1);
        sigma_n(i,j-1,16)=sigma(ind+1,1);
        sigma_n(i,j-1,13)=sigma(ind+2,1);
        sigma_n(i,j-1,14)=sigma(ind+3,1);
        sigma_n(i,j,1)=sigma(ind+2,1);
        sigma_n(i,j,2)=sigma(ind+3,1);
        sigma_n(i,j,3)=sigma(ind+4,1);
        sigma_n(i,j,4)=sigma(ind+5,1);
        sigma_n(i-1,j,5)=sigma(ind+4,1);
        sigma_n(i-1,j,6)=sigma(ind+5,1);
        sigma_n(i-1,j,7)=sigma(ind+6,1);
        sigma_n(i-1,j,8)=sigma(ind+7,1);
        sigma_n(i-1,j-1,11)=sigma(ind+6,1);
        sigma_n(i-1,j-1,12)=sigma(ind+7,1);
        ind=ind+8;
    end
    
    % East nodes 
    
    i=NN+1;
    sigma_n(i-1,j-1,9)=sigma(ind,1);
    sigma_n(i-1,j-1,10)=sigma(ind+1,1);
    sigma_n(i-1,j,5)=sigma(ind+2,1);
    sigma_n(i-1,j,6)=sigma(ind+3,1);
    sigma_n(i-1,j,7)=sigma(ind+4,1);
    sigma_n(i-1,j,8)=sigma(ind+5,1);
    sigma_n(i-1,j-1,11)=sigma(ind+4,1);
    sigma_n(i-1,j-1,12)=sigma(ind+5,1);
    ind=ind+6;
end

% NW corner node

i=1;
j=NN+1;

sigma_n(i,j-1,15)=sigma(ind,1);
sigma_n(i,j-1,16)=sigma(ind+1,1);
sigma_n(i,j-1,13)=sigma(ind+2,1);
sigma_n(i,j-1,14)=sigma(ind+3,1);
ind=ind+4;

% North nodes (j=N+1)

for i=2:NN
    sigma_n(i-1,j-1,9)=sigma(ind,1);
    sigma_n(i-1,j-1,10)=sigma(ind+1,1);
    sigma_n(i,j-1,15)=sigma(ind,1);
    sigma_n(i,j-1,16)=sigma(ind+1,1);
    sigma_n(i,j-1,13)=sigma(ind+2,1);
    sigma_n(i,j-1,14)=sigma(ind+3,1);
    sigma_n(i-1,j-1,11)=sigma(ind+4,1);
    sigma_n(i-1,j-1,12)=sigma(ind+5,1);
    ind=ind+6;
end

% NE corner node

i=NN+1;
sigma_n(i-1,j-1,9)=sigma(ind,1);
sigma_n(i-1,j-1,10)=sigma(ind+1,1);
sigma_n(i-1,j-1,11)=sigma(ind+2,1);
sigma_n(i-1,j-1,12)=sigma(ind+3,1);

return
end

% function l=side_length(i,j,k,l)
% global x y
% l=sqrt((x(i,j)-x(k,l))^2+(y(i,j)-y(k,l))^2);
% return
% end