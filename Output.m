function Output(NX,NY,phi,phi0,fid,fid_pf,cycle)
imagesc(phi(2:(NX+1),(NY+1):-1:2)'); 
colorbar
drawnow

sumdeltaphi=sum(sum(abs(phi-phi0)));
sumphi0=sum(sum(abs(phi0)));
ErrPhi=sumdeltaphi/sumphi0;
fprintf(fid,'\n%d\t%f\t%f\t%e\t%e',...
    cycle,min(phi(:)),max(phi(:)),sum(sum(phi(2:NX+1,2:NY+1))),ErrPhi);
fprintf(fid,'\r\n');

%
fprintf(fid_pf,'variables = "X", "Y","PHI" ');
fprintf(fid_pf,'\r\n');
fprintf(fid_pf,'ZONE T=''%d'',I=%d,J=%d , F=POINT',cycle,NX,NY);
fprintf(fid_pf,'\r\n');
for j=2:NY+1
    for i=2:NX+1
        fprintf(fid_pf,'\n%f\t%f\t%e ',(i-1.5)/NX,(j-1.5)/NX,phi(i,j));
        fprintf(fid_pf,'\r\n');
    end
end