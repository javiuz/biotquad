function matrix=build_matrix_sig_u

matrix=sparse(48,8);

val=1/2;
nval=-val;

% South-West corner node
rind=1:4;
cind=1:2;

b=zeros(4,2);

b(1,1)=nval;
b(3,1)=nval;
b(2,2)=nval;
b(4,2)=nval;

matrix(rind,cind)=b;

% South node 
rind=5:10;
cind=1:4;

b=zeros(6,4);
    
b(1,3)=nval;
b(2,4)=nval;
b(3,1)=val;
b(3,3)=nval;
b(4,2)=val;
b(4,4)=nval;
b(5,1)=nval;
b(6,2)=nval;

matrix(rind,cind)=b;

% South-East corner node 
rind=11:14;
cind=3:4;

b=zeros(4,2);

b(1,1)=val;
b(3,1)=nval;
b(2,2)=val;
b(4,2)=nval;

matrix(rind,cind)=b;

% West node
rind=15:20;
cind1=1:2;
cind2=5:6;

b=zeros(6,4);  

b(1,1)=nval;
b(3,1)=val;
b(2,2)=nval;
b(4,2)=val;
b(3,3)=nval;
b(5,3)=nval;
b(4,4)=nval;
b(6,4)=nval;

matrix(rind,cind1)=b(:,1:2);
matrix(rind,cind2)=b(:,3:4);
    
% Central nodes 
rind=21:28;
cind=1:8;

b=zeros(8);

b(1,1)=val;
b(7,1)=val;
b(2,2)=val;
b(8,2)=val;
b(1,3)=nval;
b(3,3)=val;
b(2,4)=nval;
b(4,4)=val;
b(5,5)=val;
b(7,5)=nval;
b(6,6)=val;
b(8,6)=nval;
b(3,7)=nval;
b(5,7)=nval;
b(4,8)=nval;
b(6,8)=nval;

matrix(rind,cind)=b;
    
% East node;
rind=29:34;
cind1=3:4;
cind2=7:8;

b=zeros(6,4);

b(1,1)=val;
b(5,1)=val;
b(2,2)=val;
b(6,2)=val;
b(3,3)=val;
b(5,3)=nval;
b(4,4)=val;
b(6,4)=nval;

matrix(rind,cind1)=b(:,1:2);
matrix(rind,cind2)=b(:,3:4);

% North-West corner node
rind=35:38;
cind=5:6;

b=zeros(4,2);

b(1,1)=nval;
b(3,1)=val;
b(2,2)=nval;
b(4,2)=val;

matrix(rind,cind)=b;

% North node
rind=39:44;
cind=5:8;

b=zeros(6,4);

b(1,1)=val;
b(5,1)=val;
b(2,2)=val;
b(6,2)=val;
b(1,3)=nval;
b(3,3)=val;
b(2,4)=nval;
b(4,4)=val;

matrix(rind,cind)=b;

% North-East corner node
rind=45:48;
cind=7:8;

b=zeros(4,2);

b(1,1)=val;
b(3,1)=val;
b(2,2)=val;
b(4,2)=val;

matrix(rind,cind)=b;
end