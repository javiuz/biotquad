function vector=build_indep_f(t)

global alpha

vector=zeros(8,1);

st=(t*alpha)/4;

vector(1)=st;
vector(3)=st;
vector(5)=st;
vector(7)=st;

end