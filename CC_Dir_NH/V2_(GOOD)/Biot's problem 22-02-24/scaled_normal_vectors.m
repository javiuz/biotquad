function [nlh,nlv]=scaled_normal_vectors(xc,yc,xh,yh,xv,yv)

if (xc-xh>0)
    nlh=[yh-yc;xc-xh];
else
    nlh=[yc-yh;xh-xc];
end

if (yc-yv>0)
    nlv=[yc-yv;xv-xc];
else
    nlv=[yv-yc;xc-xv];
end

return
end