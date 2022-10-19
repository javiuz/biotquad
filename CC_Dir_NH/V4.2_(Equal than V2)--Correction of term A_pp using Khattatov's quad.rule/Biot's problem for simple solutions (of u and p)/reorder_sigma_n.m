function [s1x,s1y,s2x,s2y]=reorder_sigma_n(sigma_n)

global NN 

N=NN;

s1x=zeros(N+1,N+1,2);
s2x=s1x;
s1y=s1x;
s2y=s1x;

% South-West corner node
i=1;
j=1;
s1y(i,j,1)=sigma_n(i,j,1);
s2y(i,j,1)=sigma_n(i,j,2);
s1x(i,j,2)=sigma_n(i,j,3);
s2x(i,j,2)=sigma_n(i,j,4);

% South nodes (j=1)
for i=2:N
    s1y(i,j,1)=sigma_n(i,j,1);
    s2y(i,j,1)=sigma_n(i,j,2);
    s1x(i,j,2)=sigma_n(i,j,3);
    s2x(i,j,2)=sigma_n(i,j,4);
    s1y(i,j,2)=sigma_n(i-1,j,7);
    s2y(i,j,2)=sigma_n(i-1,j,8);
end

% South-East corner node (j=1)
i=N+1;
s1x(i,j,2)=sigma_n(i-1,j,5);
s2x(i,j,2)=sigma_n(i-1,j,6);
s1y(i,j,2)=sigma_n(i-1,j,7);
s2y(i,j,2)=sigma_n(i-1,j,8);

for j=2:N

    % West nodes
    i=1;
    s1x(i,j,1)=sigma_n(i,j-1,15);
    s2x(i,j,1)=sigma_n(i,j-1,16);
    s1y(i,j,1)=sigma_n(i,j,1);
    s2y(i,j,1)=sigma_n(i,j,2);
    s1x(i,j,2)=sigma_n(i,j,3);
    s2x(i,j,2)=sigma_n(i,j,4);
    
    % Central nodes 
    for i=2:N
        s1x(i,j,1)=sigma_n(i,j-1,15);
        s2x(i,j,1)=sigma_n(i,j-1,16);
        s1y(i,j,1)=sigma_n(i,j,1);
        s2y(i,j,1)=sigma_n(i,j,2);
        s1x(i,j,2)=sigma_n(i,j,3);
        s2x(i,j,2)=sigma_n(i,j,4);
        s1y(i,j,2)=sigma_n(i-1,j,7);
        s2y(i,j,2)=sigma_n(i-1,j,8);
    end
    
    % East nodes
    i=N+1;
    s1x(i,j,1)=sigma_n(i-1,j-1,9);
    s2x(i,j,1)=sigma_n(i-1,j-1,10);
    s1x(i,j,2)=sigma_n(i-1,j,5);
    s2x(i,j,2)=sigma_n(i-1,j,6);
    s1y(i,j,2)=sigma_n(i-1,j,7);
    s2y(i,j,2)=sigma_n(i-1,j,8);
end

% North-West corner node
i=1;
j=N+1;
s1x(i,j,1)=sigma_n(i,j-1,15);
s2x(i,j,1)=sigma_n(i,j-1,16);
s1y(i,j,1)=sigma_n(i,j-1,13);
s2y(i,j,1)=sigma_n(i,j-1,14);

% North nodes (j=N+1)
for i=2:N
    s1x(i,j,1)=sigma_n(i,j-1,15);
    s2x(i,j,1)=sigma_n(i,j-1,16);
    s1y(i,j,1)=sigma_n(i,j-1,13);
    s2y(i,j,1)=sigma_n(i,j-1,14);
    s1y(i,j,2)=sigma_n(i-1,j-1,11);
    s2y(i,j,2)=sigma_n(i-1,j-1,12);
end

% North-East corner node (j=N+1)
i=N+1;
s1x(i,j,1)=sigma_n(i-1,j-1,9);
s2x(i,j,1)=sigma_n(i-1,j-1,10);
s1y(i,j,2)=sigma_n(i-1,j-1,11);
s2y(i,j,2)=sigma_n(i-1,j-1,12);

end