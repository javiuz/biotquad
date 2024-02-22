function matrix=build_matrix_p_tau_manual

global alpha lambda mu

coef=alpha/(16*(lambda+mu));

matrix=sparse(48,4);

% South-West corner node
rind=1:4;
cind=1;

d=zeros(4,1);

d(2,1)=1;
d(3,1)=1;

matrix(rind,cind)=d;

% South node 
rind=5:10;
cind=1:2;

d=zeros(6,2);
    
d(3,1)=1;
d(6,1)=1;
d(2,2)=1;
d(3,2)=1;

matrix(rind,cind)=d;

% South-East corner node 
rind=11:14;
cind=2;

d=zeros(4,1);

d(1,1)=1;
d(4,1)=1;

matrix(rind,cind)=d;

% West node
rind=15:20;
cind1=1;
cind2=3;

d=zeros(6,2);  

d(1,1)=1;
d(4,1)=1;
d(4,2)=1;
d(5,2)=1;

matrix(rind,cind1)=d(:,1);
matrix(rind,cind2)=d(:,2);
    
% Central nodes 
rind=21:28;
cind=1:4;

d=zeros(8,4);

d(1,1)=1;
d(8,1)=1;
d(1,2)=1;
d(4,2)=1;
d(5,3)=1;
d(8,3)=1;
d(4,4)=1;
d(5,4)=1;

matrix(rind,cind)=d;
    
% East node;
rind=29:34;
cind1=2;
cind2=4;

d=zeros(6,2);

d(1,1)=1;
d(6,1)=1;
d(3,2)=1;
d(6,2)=1;

matrix(rind,cind1)=d(:,1);
matrix(rind,cind2)=d(:,2);

% North-West corner node
rind=35:38;
cind=3;

d=zeros(4,1);

d(1,1)=1;
d(4,1)=1;

matrix(rind,cind)=d;

% North node
rind=39:44;
cind=3:4;

d=zeros(6,2);

d(1,1)=1;
d(6,1)=1;
d(1,2)=1;
d(4,2)=1;

matrix(rind,cind)=d;

% North-East corner node
rind=45:48;
cind=4;

d=zeros(4,1);

d(1,1)=1;
d(4,1)=1;

matrix(rind,cind)=d;
matrix=coef*matrix;
end