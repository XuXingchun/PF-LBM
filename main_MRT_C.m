format long
% Lattice parameters
NX=500;NY=500;R=100;
X0=NX/2+1;Y0=NY/2+1;
dx=1.0;dy=1.0;dt=dx;
Lx=dx*NX;

e=[0 0;1 0;0 1;-1 0;0 -1;1 1;-1 1;-1 -1;1 -1]; % D2Q9
w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]; % weight
M=[1 1 1 1 1 1 1 1 1; -4 -1 -1 -1 -1 2 2 2 2; 4 -2 -2 -2 -2 1 1 1 1;...
    0 1 0 -1 0 1 -1 -1 1; 0 -2 0 2 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1;...
    0 0 -2 0 2 1 1 -1 -1; 0 1 -1 1 -1 0 0 0 0; 0 0 0 0 0 1 -1 1 -1];
%
phih=1;phil=-1;
Wc=2;
U0=0.25;
Tf=1.25* NX /U0;
maxT  =2* Tf;
% -------relaxation matrix-------%-----------------------------------------
Mphi=0.1;
se=0.6;
tau_c=3 * Mphi  + 0.5;
S= [1.0 se se 1/tau_c 1/tau_c 1/tau_c 1/tau_c 1.0 1.0];    
S=diag(S);
inv_M_S = M\S*M;
%--------------------------------------------------------------------------
filename=sprintf('MRTC_His_U0%1.3f_M%1.4f_se%1.2f.txt',U0,Mphi,se);
fid=fopen(filename,'w+');
fprintf(fid,'t  phi_min  phi_max  Mass_C  ErrPhi \r\n');

file_pf=sprintf('MRTC_Res_U0%1.3f_M%1.4f_se%1.2f.dat',U0,Mphi,se);
fid_pf=fopen(file_pf,'w');

%% initialize
ux=zeros(NX+2,NY+2);
uy=zeros(NX+2,NY+2);
phi=zeros(NX+2,NY+2); 

%
for j=1:NY+2
    for i=1:NX+2
        Ri=sqrt((i-(X0+0.5))^2+(j-(Y0+0.5))^2);
        phi(i,j)=(phih+phil)/2-(phih-phil)/2*tanh(2*(R-Ri)/Wc);
        
        ux(i,j) = - U0 *sin(4* pi * ((i - 1.5)/Lx +0.5)) *sin(4* pi * ((j - 1.5)/Lx +0.5));
        uy(i,j) = - U0 *cos(4* pi * ((i - 1.5)/Lx +0.5)) *cos(4* pi * ((j - 1.5)/Lx +0.5));
    end
end
phi0=phi;ux0=ux;uy0=uy;

[heq, hlpC]=Calculate_heq_G_MRT_C(NX,NY,Wc, e, w, ux, uy, phi, phil, phih);
h=heq;

Output(NX,NY,phi,phi0,fid,fid_pf,0)
%% Ë²Ì¬·ÖÎö
for cycle = 1:maxT
    % collsion
    hneq=heq-h;
    hneq=reshape(hneq,(NX+2)*(NY+2),9);
    omehneq=(inv_M_S*hneq')';
    omehneq=reshape(omehneq,(NX+2),(NY+2),9);
    h=h + omehneq + hlpC;
    
    %periodic B.C.
    h(1   ,:,:) = h(NX+1,:,:);	
    h(NX+2,:,:) = h(2   ,:,:);	
    h(:, 1  ,:) = h(:,NY+1,:);	
    h(:,NY+2,:) = h(:,2   ,:);

    %Streaming
    htmp=h;
    for k=1:9
        h(:,:,k)=circshift(htmp(:,:,k),[e(k,1),e(k,2)]); %
    end
    
    %Macroscopic_Properties
    phi = sum( h(:,:,:), 3);
    
    phi(:,  1 ) = phi(:,NY+1  );
    phi(:,NY+2) = phi(:, 2  );
    phi(  1 ,:) = phi(NX+1,:);
    phi(NX+2,:) = phi( 2  ,:);
    %
    ux = ux0.*cos(pi*cycle/Tf);
    uy = uy0.*cos(pi*cycle/Tf);
    
    [heq, hlpC]=Calculate_heq_G_MRT_C(NX,NY,Wc, e, w, ux, uy, phi, phil, phih);
    %
    if mod( cycle,fix(Tf/4))==0
        cycle
        max(phi(:))   
        Output(NX,NY,phi,phi0,fid,fid_pf,cycle)
    end
end