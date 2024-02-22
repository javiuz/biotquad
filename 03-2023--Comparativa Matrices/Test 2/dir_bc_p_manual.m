function vector=dir_bc_p_manual(t)

vector=zeros(24,1);

% South-West corner node
vector(1)=t/8;

% South node 
vector(3)=3*t/8;
vector(5)=vector(1);

% South-East corner node 
vector(6)=-t/2;
vector(7)=3*t/8;

% West node: =0
    
% Central nodes : No hay
    
% East node;
vector(15)=-t/2;
vector(16)=-t/2;

% North-West corner node
vector(19)=-t/8;

% North node
vector(21)=-3*t/8;
vector(22)=-t/8;

% North-East corner node
vector(23)=-t/2;
vector(24)=-3*t/8;

end