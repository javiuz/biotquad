function matrix=build_matrix_sig_gamma_manual

matrix=sparse(48,9);

% South-West corner node
rind=1:4;
cind=1;

c=zeros(4,1);

c(1,1)=1/8;
c(4,1)=-1/8;

matrix(rind,cind)=c;

% South node 
rind=5:10;
cind=2;

c=zeros(6,1);
    
c(1,1)=1/8;
c(4,1)=-1/4;
c(5,1)=1/8;

matrix(rind,cind)=c;

% South-East corner node 
rind=11:14;
cind=3;

c=zeros(4,1);

c(2,1)=-1/8;
c(3,1)=1/8;

matrix(rind,cind)=c;

% West node
rind=15:20;
cind=4;

c=zeros(6,1);  

c(2,1)=-1/8;
c(3,1)=1/4;
c(6,1)=-1/8;

matrix(rind,cind)=c;
    
% Central nodes 
rind=21:28;
cind=5;

c=zeros(8,1);

c(2,1)=-1/4;
c(3,1)=1/4;
c(6,1)=-1/4;
c(7,1)=1/4;

matrix(rind,cind)=c;
    
% East node;
rind=29:34;
cind=6;

c=zeros(6,1);

c(2,1)=-1/8;
c(4,1)=-1/8;
c(5,1)=1/4;

matrix(rind,cind)=c;

% North-West corner node
rind=35:38;
cind=7;

c=zeros(4,1);

c(2,1)=-1/8;
c(3,1)=1/8;

matrix(rind,cind)=c;

% North node
rind=39:44;
cind=8;

c=zeros(6,1);

c(2,1)=-1/4;
c(3,1)=1/8;
c(5,1)=1/8;

matrix(rind,cind)=c;

% North-East corner node
rind=45:48;
cind=9;

c=zeros(4,1);

c(2,1)=-1/8;
c(3,1)=1/8;

matrix(rind,cind)=c;
end