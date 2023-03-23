function matrix=build_matrix_z_z

matrix=sparse(24,24);

% South-West corner node
ind=1:2;

w=zeros(2);

w(1,1)=1/4;
w(2,2)=1/4;

matrix(ind,ind)=w;

% South node 
ind=3:5;

w=zeros(3);
    
w(1,1)=1/4;
w(2,2)=1/2;
w(3,3)=1/4;

matrix(ind,ind)=w;

% South-East corner node 
ind=6:7;

w=zeros(2);

w(1,1)=1/4;
w(2,2)=1/4;

matrix(ind,ind)=w;

% West node
ind=8:10;

w=zeros(3);  

w(1,1)=1/4;
w(2,2)=1/2;
w(3,3)=1/4;    

matrix(ind,ind)=w;
    
% Central nodes 
ind=11:14;

w=zeros(4);

w(1,1)=1/2;
w(2,2)=w(1,1);
w(3,3)=w(1,1);
w(4,4)=w(1,1);

matrix(ind,ind)=w;
    
% East node
ind=15:17;

w=zeros(3);

w(1,1)=1/4;
w(2,2)=1/4;
w(3,3)=1/2;

matrix(ind,ind)=w;

% North-West corner node
ind=18:19;

w=zeros(2);

w(1,1)=1/4;
w(2,2)=1/4;

matrix(ind,ind)=w;

% North node
ind=20:22;

w=zeros(3);

w(1,1)=1/2;
w(2,2)=1/4;
w(3,3)=1/4;

matrix(ind,ind)=w;

% North-East corner node
ind=23:24;

w=zeros(2);

w(1,1)=1/4;
w(2,2)=1/4;

matrix(ind,ind)=w;
end