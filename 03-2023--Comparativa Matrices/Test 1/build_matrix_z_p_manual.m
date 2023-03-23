function matrix=build_matrix_z_p_manual

matrix=sparse(24,4);

val=1/2;
nval=-val;

% South-West corner node
rind=1:2;
cind=1;

f=zeros(2,1);

f(1,1)=val;
f(2,1)=val;

matrix(rind,cind)=f;

% South node 
rind=3:5;
cind=1:2;

f=zeros(3,2);
    
f(2,1)=nval;
f(3,1)=val;
f(1,2)=val;
f(2,2)=val;

matrix(rind,cind)=f;

% South-East corner node 
rind=6:7;
cind=2;

f=zeros(2,1);

f(1,1)=nval;
f(2,1)=val;

matrix(rind,cind)=f;

% West node
rind=8:10;
cind1=1;
cind2=3;

f=zeros(3,2);  

f(1,1)=val;
f(2,1)=val;
f(2,2)=nval;
f(3,2)=val;

matrix(rind,cind1)=f(:,1);
matrix(rind,cind2)=f(:,2);
    
% Central nodes 
rind=11:14;
cind=1:4;

f=zeros(4);

f(1,1)=val;
f(4,1)=val;
f(1,2)=nval;
f(2,2)=val;
f(3,3)=val;
f(4,3)=nval;
f(2,4)=nval;
f(3,4)=nval;

matrix(rind,cind)=f;
    
% East node;
rind=15:17;
cind1=2;
cind2=4;

f=zeros(3,2);

f(1,1)=nval;
f(3,1)=val;
f(2,2)=nval;
f(3,2)=nval;

matrix(rind,cind1)=f(:,1);
matrix(rind,cind2)=f(:,2);

% North-West corner node
rind=18:19;
cind=3;

f=zeros(2,1);

f(1,1)=val;
f(2,1)=nval;

matrix(rind,cind)=f;

% North node
rind=20:22;
cind=3:4;

f=zeros(3,2);

f(1,1)=nval;
f(3,1)=nval;
f(1,2)=val;
f(2,2)=nval;

matrix(rind,cind)=f;

% North-East corner node
rind=23:24;
cind=4;

f=zeros(2,1);

f(1,1)=nval;
f(2,1)=nval;

matrix(rind,cind)=f;
end