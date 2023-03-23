function vector=dir_bc_u

vector=zeros(48,1);

% South-West corner node
vector(1)=-1/8;
vector(4)=vector(1);

% South node 
vector(5)=-3/8;
vector(9)=vector(1);

% South-East corner node 
vector(11)=1/2;
vector(12)=1/8;
vector(13)=-3/8;

% West node
vector(16)=-1/8;
vector(20)=-3/8;
    
% Central nodes : No hay
    
% East node;
vector(29)=1/2;
vector(30)=1/8;
vector(31)=1/2;
vector(32)=3/8;

% North-West corner node
vector(36)=-3/8;
vector(37)=1/8;
vector(38)=1/2;

% North node
vector(41)=3/8;
vector(42)=1/2;
vector(43)=1/8;
vector(44)=1/2;

% North-East corner node
vector(45)=1/2;
vector(46)=3/8;
vector(47)=3/8;
vector(48)=1/2;

end