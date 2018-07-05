load w_snap_sample.mat
load FTLE_snap_new.mat

nt=size(w_snap,2);
ny=513;nx=257;

dx=1/256;
dy=dx;
W=simpson_2D_weights(nx,ny,dx,dy);
alpha=1;

for i=1:nt
    for j=i:nt
        w1=reshape(w_snap(:,i),[513 257]);
        w2=reshape(w_snap(:,j),[513 257]);
        
        l1=reshape(ftle_snap(:,i),[513 257]);
        l2=reshape(ftle_snap(:,j),[513 257]);
        C(i,j)=sum(sum(W.*(w1.*w2 + alpha*l1.*l2)));
        %C(i,j)=sum(sum(W.*(w1.*w2)));
        C(j,i)=C(i,j);
    end
    i,
end
C=C./nt;
% C=w_snap'*M*w_snap;
% C=C./nt;

[V,D]=eig(C);
[D,I]=sort(diag(D),'descend');
V=V(:,I);

nModes=100;
N=ny*nx;
wPOD=zeros(N,nModes);
ftlePOD=zeros(N,nModes);

for i=1:nModes
    for j=1:nt
        %wPOD(:,:,i)=wPOD(:,:,i)+V(j,i)*w(:,:,j);
        wPOD(:,i)   = wPOD(:,i) + V(j,i)*w_snap(:,j);
        %ftlePOD(:,i)=ftlePOD(:,i)+V(j,i)*ftle_snap(:,j);
    end
    %dd(j)=sum(V(:,j).*V(:,j));
    wPOD(:,i)=wPOD(:,i)./sqrt(nt*D(i));
    %ftlePOD(:,i)=ftlePOD(:,i)./sqrt(nt*D(i));
    %----Normalize
    %modeFactor = 1./sqrt(nt*D(i));
    %wPOD(:,:,i)=wPOD(:,:,i)*modeFactor;
end



Id=zeros(nModes,nModes);
for i=1:nModes
    for j=1:nModes
        w1=reshape(wPOD(:,i),[513 257]);
        w2=reshape(wPOD(:,j),[513 257]);
        l1=reshape(ftlePOD(:,i),[513 257]);
        l2=reshape(ftlePOD(:,j),[513 257]);
        Id(i,j)=sum(sum(W.*(w1.*w2+alpha*l1.*l2)));
        %Id(i,j)=sum(sum(W.*(w1.*w2)));
        %Id(i,j)=sum(sum(W.*(wPOD(:,:,i).*wPOD(:,:,j))));
    end
end
cond(Id),

load FE_matrix_top.mat
 elem=Mesh_trig_node;
 nodes=Mesh_node_xy;
 FEM.elem=elem;
 FEM.nodes=nodes;
 
 psi_d_name='d_func';
 xmin=0;xmax=1;
 ymin=0;ymax=2;
  
 domain = [xmin,xmax,ymin,ymax];
 bc_index_q=[0 0 0 0];
 
 GDOF_q.P_g=FEM.nodes';
 GDOF_q.T_g=FEM.elem';
 [node_type_q,DOF_q] = global_dof_p_2D_fe_Lagrange(domain,bc_index_q,GDOF_q);
 
 q_gd = FE_u_gd_2D_Lagrange(psi_d_name,node_type_q,GDOF_q);


tM=M(DOF_q,DOF_q);
S=M_xx+M_yy;
tS=S(DOF_q,DOF_q);

[N, k]=size(wPOD);

POD_psi=zeros(N,nModes);

for i=1:nModes    
    vor=wPOD(:,i);
    tvor=vor(DOF_q);
    tpsi=tS\(tM*tvor);
    psi=FE_tu_fe_to_u_fe_2D_Lagrange(tpsi,q_gd,node_type_q,DOF_q);    
    POD_psi(:,i)=psi;
end








