function [heq, hlpC]=Calculate_heq_G_MRT_C(NX,NY,Wc, e, w, ux, uy, phi, phil, phih)
DcDx=zeros(NX+2,NY+2);
DcDy=zeros(NX+2,NY+2);
Gamma=zeros(NX+2,NY+2,9);
heq=zeros(NX+2,NY+2,9);
hlpC=zeros(NX+2,NY+2,9);

for j = 2: NY+1
    for i = 2: NX+1
        DcDx(i,j) = (phi(i+1,j  ) - phi(i-1,j  ))/3 ...
            + ( phi(i+1,j-1) + phi(i+1,j+1) - phi(i-1,j-1) - phi(i-1,j+1))/12;
        DcDy(i,j) = (phi(i  ,j+1) - phi(i  ,j-1))/3 ...
            + ( phi(i-1,j+1) + phi(i+1,j+1) - phi(i-1,j-1) - phi(i+1,j-1))/12;  
    end
end
tmp = ( DcDx.^2 + DcDy.^2 + 1e-32 ).^0.5;
ni = DcDx ./ tmp;
nj = DcDy ./ tmp;
%--------------------------------------------------------------------------------
uv=ux.*ux+uy.*uy;
for k=1:9
    eF = 4.*(phi-phih).*(phi -phil)./Wc/(phil-phih).*(e(k,1).*ni + e(k,2).*nj);
    hlpC(:,:,k) = w(k) .* eF;
    eu=e(k,1).*ux+e(k,2).*uy;
    Gamma(:,:,k)=w(k).*(1+3.*eu+4.5* eu.*eu-1.5.*uv);
    heq(:,:,k) = phi.*Gamma(:,:,k) - 0.5 * hlpC(:,:,k);
end
%--------------------------------------------------------------------------------