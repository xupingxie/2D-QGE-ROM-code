
%function [POD_w,POD_psi,Diag_S]=compute_POD_w_psi_basis(u_snap,M,S,r,alpha,FEM)
%load FE_matrix.mat
%load vorticity_flucuation_snap.mat

%% ---compute POD_vor by snapshots
%u_snap=vorprime_snap(:,101:801);
alpha=0;

[POD_vor,Diag_S,d,CumEng,CumEng_ratio,rows,dim] = PODbasis_QGE(u_snap,M);
%if alpha==0
%    POD_w = POD_vor(:,1:r);
%else
    POD_w = POD_vor(rows(1):rows(2)-1,1:r);
%end

%mr= cond(POD_w(:,1:(end-1))'*M*POD_w(:,1:(end-1)));

%if mr==1
%    fprintf('It is correct\n')
%else
%    fprintf('WRONG!\n')
%end

%% ----compute POD_psi by Possion assumption


psi_d_name='d_func';
 xmin=0;xmax=1;
 ymin=0;ymax=2;
 
domain = [xmin,xmax,ymin,ymax];
bc_index_q=[0 0 0 0];

GDOF_q.P_g=FEM.nodes';
GDOF_q.T_g=FEM.elem';
[node_type_q,DOF_q] = global_dof_p_2D_fe_Lagrange(domain,bc_index_q,GDOF_q);

q_gd = FE_u_gd_2D_Lagrange(psi_d_name,node_type_q,GDOF_q);
%FE_tu_fe_to_u_fe_2D_Lagrange(tp0,q_gd,node_type_q,DOF_q);
tM=M(DOF_q,DOF_q);
%S=M_xx+M_yy;
tS=S(DOF_q,DOF_q);

[N, k]=size(POD_w);

POD_psi=zeros(N,r);

for i=1:r    
    vor=POD_w(:,i);
    tvor=vor(DOF_q);
    tpsi=tS\(tM*tvor);
    psi=FE_tu_fe_to_u_fe_2D_Lagrange(tpsi,q_gd,node_type_q,DOF_q);    
    POD_psi(:,i)=psi;
end
