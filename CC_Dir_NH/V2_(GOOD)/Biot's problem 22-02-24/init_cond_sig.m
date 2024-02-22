function sig=init_cond_sig

global x y NN

N=NN;

sig=zeros(8*N*(N+1),1);

% SW corner node
i=1;
j=1;

v1=1;
v2=2;
h1=3;
h2=4;

xc=x(i,j);
yc=y(i,j);
xh=x(i+1,j);
yh=y(i+1,j);
xv=x(i,j+1);
yv=y(i,j+1);

[nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

e=exact_sol_sigma(xc,yc,0);
sig(h1:h2,1)=e*nlv;
sig(v1:v2,1)=e*nlh;

% S nodes
% j=1;
for i=2:N

     % Clculos sobre el elemento W
    v1=9+(i-2)*6;
    v2=v1+1;
    h1=7+(i-2)*6;
    h2=h1+1;

    xc=x(i,j);
    yc=y(i,j);
    xh=x(i-1,j);
    yh=y(i-1,j);
    xv=x(i,j+1);
    yv=y(i,j+1);

    [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

    e=exact_sol_sigma(xc,yc,0);
    sig(h1:h2,1)=e*nlv;
    sig(v1:v2,1)=e*nlh;

    % Clculos sobre el elemento E
    v1=5+(i-2)*6;
    v2=v1+1;

    % xc,yc,xv,yv no cambian
    xh=x(i+1,j);
    yh=y(i+1,j);
    
[nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

e=exact_sol_sigma(xc,yc,0);
%sig(h1:h2,1)=e*nlv;
sig(v1:v2,1)=e*nlh;

end

% SE corner node
% j=1;
i=N+1;

h1=5+(N-1)*6;
h2=h1+1;
v1=h1+2;
v2=h1+3;

xc=x(i,j);
yc=y(i,j);
xh=x(i-1,j);
yh=y(i-1,j);
xv=x(i,j+1);
yv=y(i,j+1);

[nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

e=exact_sol_sigma(xc,yc,0);
sig(h1:h2,1)=e*nlv;
sig(v1:v2,1)=e*nlh;

for j=2:N
    % West nodes:
    i=1;
    
    % Clculos sobre elemento S

    h1=6*N+3+(j-2)*(8*N+4);
    h2=h1+1;
    v1=h1+2;
    v2=h1+3;

    xc=x(i,j);
    yc=y(i,j);
    xh=x(i+1,j);
    yh=y(i+1,j);
    xv=x(i,j-1);
    yv=y(i,j-1);
    
    [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

    e=exact_sol_sigma(xc,yc,0);
    sig(h1:h2,1)=e*nlv;
    sig(v1:v2,1)=e*nlh;

    % Clculos sobre elemento N

    v1=6*N+5+(j-2)*(8*N+4);
    v2=v1+1;
    h1=v1+2;
    h2=v1+3;

    % xc,yc,xh,yh no cambian
    xv=x(i,j+1);
    yv=y(i,j+1);

    [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

    e=exact_sol_sigma(xc,yc,0);
    sig(h1:h2,1)=e*nlv;
    %sig(v1:v2,1)=e*nlh;

    % Central nodes:
    for i=2:N
       
        % Clculos sobre elemento SW

        h1=6*N+9+(j-2)*(8*N+4)+(i-2)*8;
        h2=h1+1;
        v1=h1+6;
        v2=h1+7;

        xc=x(i,j);
        yc=y(i,j);
        xh=x(i-1,j);
        yh=y(i-1,j);
        xv=x(i,j-1);
        yv=y(i,j-1);

        [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

        e=exact_sol_sigma(xc,yc,0);
        sig(h1:h2,1)=e*nlv;
        sig(v1:v2,1)=e*nlh;

        % Clculos sobre elemento SE
        % h1, h2 no cambian
        v1=h1+2;
        v2=h1+3;

        % xc, yc, xv, yv no cambian
        xh=x(i+1,j);
        yh=y(i+1,j);
        
        [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

        e=exact_sol_sigma(xc,yc,0);
        %sig(h1:h2,1)=e*nlv;
        sig(v1:v2,1)=e*nlh;

        % Clculos sobre el elemento NE
        % v1, v2 no cambian
        h1=v2+1;
        h2=v2+2;
        
        % xc, yc, xh, yh no cambian
        xv=x(i,j+1);
        yv=y(i,j+1);

        [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

        e=exact_sol_sigma(xc,yc,0);
        sig(h1:h2,1)=e*nlv;
        %sig(v1:v2,1)=e*nlh;

        % Sobre el elemento NW no calculamos nada, ya están todos
        % calculados
    end

    % East nodes
    i=N+1;

    % Clculos sobre elemento S

    h1=6*N+3+(j-1)*(8*N+4)-6;
    h2=h1+1;
    v1=h1+4;
    v2=h1+5;

    xc=x(i,j);
    yc=y(i,j);
    xh=x(i-1,j);
    yh=y(i-1,j);
    xv=x(i,j-1);
    yv=y(i,j-1);

    [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

    e=exact_sol_sigma(xc,yc,0);
    sig(h1:h2,1)=e*nlv;
    sig(v1:v2,1)=e*nlh;

    % Clculos sobre elemento N
    % v1 y v2 no cambian
    h1=v1-2;
    h2=v1-1;

    % xc, yc, xh, yh no cambian
    xv=x(i,j+1);
    yv=y(i,j+1);

    [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

    e=exact_sol_sigma(xc,yc,0);
    sig(h1:h2,1)=e*nlv;
    %sig(v1:v2,1)=e*nlh;
end

% NW corner node
i=1;
j=N+1;

h1=6*N+2+(N-1)*(8*N+4)+1;
h2=h1+1;
v1=h1+2;
v2=h1+3;

xc=x(i,j);
yc=y(i,j);
xh=x(i+1,j);
yh=y(i+1,j);
xv=x(i,j-1);
yv=y(i,j-1);

[nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

e=exact_sol_sigma(xc,yc,0);
sig(h1:h2,1)=e*nlv;
sig(v1:v2,1)=e*nlh;

% N nodes
% j=N+1;
for i=2:N
   
    % Clculos sobre el elemento W

    h1=6*N+2+(N-1)*(8*N+4)+5+(i-2)*6;
    h2=h1+1;
    v1=h1+4;
    v2=h1+5;

    xc=x(i,j);
    yc=y(i,j);
    xh=x(i-1,j);
    yh=y(i-1,j);
    xv=x(i,j-1);
    yv=y(i,j-1);

    [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

    e=exact_sol_sigma(xc,yc,0);
    sig(h1:h2,1)=e*nlv;
    sig(v1:v2,1)=e*nlh;

    % Clculos sobre el elemento E
    % h1, h2 no cambian
    v1=h1+2;
    v2=h1+3;

    % xc,yc,xv,yv no cambian
    xh=x(i+1,j);
    yh=y(i+1,j);
    
    [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

    e=exact_sol_sigma(xc,yc,0);
    %sig(h1:h2,1)=e*nlv;
    sig(v1:v2,1)=e*nlh;
end

% NE corner node
% j=N+1;
i=N+1;

h1=8*N*(N+1)-3;
h2=h1+1;
v1=h1+2;
v2=h1+3;

xc=x(i,j);
yc=y(i,j);
xh=x(i-1,j);
yh=y(i-1,j);
xv=x(i,j-1);
yv=y(i,j-1);

[nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv);

e=exact_sol_sigma(xc,yc,0);
sig(h1:h2,1)=e*nlv;
sig(v1:v2,1)=e*nlh;

return
end