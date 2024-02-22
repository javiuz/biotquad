function matrix=build_matrix_sig_sig_manual

global lambda mu

matrix=sparse(48,48);

d16mlm=16*mu*(lambda + mu);
d8mlm=d16mlm/2;
d8m=8*mu;
d4m=d8m/2;

l2m=lambda + 2*mu;

% South-West corner node
ind=1:4;

a=zeros(4);

a(1,1)=1/d8m;
a(2,2)=l2m/d16mlm;
a(2,3)=-lambda/d16mlm;
a(3,2)=a(2,3);
a(3,3)=a(2,2);
a(4,4)=a(1,1);

matrix(ind,ind)=a;

% South node 
ind=5:10;

a=zeros(6);
    
a(1,1)=1/d8m;
a(2,2)=l2m/d16mlm;
a(2,3)=-lambda/d16mlm;
a(3,2)=a(2,3);
a(3,3)=l2m/d8mlm;
a(3,6)=a(2,3);
a(4,4)=1/d4m;
a(5,5)=a(1,1);
a(6,3)=a(3,6);
a(6,6)=a(2,2);

matrix(ind,ind)=a;

% South-East corner node 
ind=11:14;

a=zeros(4);

a(1,1)=l2m/d16mlm;
a(1,4)=-lambda /d16mlm;
a(2,2)=1/d8m;
a(3,3)=a(2,2);
a(4,1)=a(1,4);
a(4,4)=a(1,1);

matrix(ind,ind)=a;

% West node
ind=15:20;

a=zeros(6);  

a(1,1)=l2m/d16mlm;
a(1,4)=-lambda/d16mlm;
a(2,2)=1/d8m;
a(3,3)=1/d4m;    
a(4,1)=a(1,4);
a(4,4)=l2m/d8mlm;
a(4,5)=a(1,4);
a(5,4)=a(4,5);
a(5,5)=a(1,1);
a(6,6)=a(2,2);

matrix(ind,ind)=a;
    
% Central nodes 
ind=21:28;

a=zeros(8);

a(1,1)=l2m/d8mlm;
a(1,4)=-lambda/d16mlm;
a(1,8)=a(1,4);
a(2,2)=1/d4m;
a(3,3)=a(2,2);
a(4,1)=a(1,4);
a(4,4)=a(1,1);
a(4,5)=a(1,8);
a(5,4)=a(4,5);
a(5,5)=a(1,1);
a(5,8)=a(1,8);
a(6,6)=a(2,2);
a(7,7)=a(2,2);
a(8,1)=a(1,8);
a(8,5)=a(5,8);
a(8,8)=a(1,1);

matrix(ind,ind)=a;
    
% East node;
ind=29:34;

a=zeros(6);

a(1,1)=l2m/d16mlm;
a(1,6)=-lambda/d16mlm;
a(2,2)=1/d8m;
a(3,3)=a(1,1);
a(3,6)=a(1,6);
a(4,4)=a(2,2);
a(5,5)=1/d4m;
a(6,1)=a(1,6);
a(6,3)=a(3,6);
a(6,6)=l2m/d8mlm;

matrix(ind,ind)=a;

% North-West corner node
ind=35:38;

a=zeros(4);

a(1,1)=l2m/d16mlm;
a(1,4)=-lambda/d16mlm;
a(2,2)=1/d8m;
a(3,3)=a(2,2);
a(4,1)=a(1,4);
a(4,4)=a(1,1);

matrix(ind,ind)=a;

% North node
ind=39:44;

a=zeros(6);

a(1,1)=l2m/d8mlm;
a(1,4)=-lambda/d16mlm;
a(1,6)=a(1,4);
a(2,2)=1/d4m;
a(3,3)=1/d8m;
a(4,1)= a(1,4);
a(4,4)=l2m/d16mlm;
a(5,5)=a(3,3);
a(6,1)=a(1,6);
a(6,6)=a(4,4);

matrix(ind,ind)=a;

% North-East corner node
ind=45:48;

a=zeros(4);

a(1,1)=l2m/d16mlm;
a(1,4)=-lambda/d16mlm;
a(2,2)=1/d8m;
a(3,3)=a(2,2);
a(4,1)=a(1,4);
a(4,4)=a(1,1);

matrix(ind,ind)=a;
end