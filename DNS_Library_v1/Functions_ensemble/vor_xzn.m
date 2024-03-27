function [gx,gz] = vor_xzn(u0,v0,w0,Fn)
% x- and z- normal Vorticity fields 

% global a b N NX MZ D1x D2x D1z D2z DY D2 DYF D2F yE xE zE dy dx dz A B L 

dydw=difY_Fn(w0,1,Fn);
dzdv=difZ_Fn(v0,1,Fn);

dxdv=difX_Fn(v0,1,Fn);
dydu=difY_Fn(u0,1,Fn);

gx=dydw-dzdv;
gz=dxdv-dydu;


end

