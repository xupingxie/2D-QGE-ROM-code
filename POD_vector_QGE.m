
nt=size(w_snap,2);
ny=513;nx=257;


dx=1/256;
dy=dx;
W=simpson_2D_weights(nx,ny,dx,dy);

w=zeros(513,257,nt);
C=zeros(nt,nt);
for k=1:nt
    w(:,:,k)=reshape(w_snap(:,k),[513 257]);
end
for i=1:nt
    for j=i:nt
        %w(:,:,i)=reshape(w_snap(:,i),[513 257]);
        %w(:,:,j)=reshape(w_snap(:,j),[513 257]);
        C(i,j)=sum(sum(W.*(w(:,:,i).*w(:,:,j))));
        C(j,i)=C(i,j);
    end
end
C=C./nt;

[V,D]=eig(C);
[D,I]=sort(diag(D),'descend');
V=V(:,I);

nModes=30;
wPOD=zeros(ny,nx,nModes);
for i=1:nModes
    for j=1:nt
        wPOD(:,:,i)=wPOD(:,:,i)+V(j,i)*w(:,:,j);
    end
    
    %----Normalize
    modeFactor = 1./sqrt(nt*D(i));
    wPOD(:,:,i)=wPOD(:,:,i)*modeFactor;
end


Id=zeros(nModes,nModes);
for i=1:nModes
    for j=1:nModes
        Id(i,j)=sum(sum(W.*(wPOD(:,:,i).*wPOD(:,:,j))));
    end
end