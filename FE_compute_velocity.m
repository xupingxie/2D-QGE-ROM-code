
% skip=10;
% idx=1:skip:257;
% idy=1:skip:513;
% k=1001;
% p=reshape(psi_snap(:,k),[513 257]);
% u=reshape(u_test(:,k),[513 257]);
% v=reshape(v_test(:,k),[513 257]);
% hold on
% contour(xs,ys,p,20);
% quiver(xs(idx),ys(idy),u(idy,idx),v(idy,idx));
% hold off

psi_d_name='d_func';
xmin=0;xmax=1;
ymin=0;ymax=2;
domain = [xmin,xmax,ymin,ymax];
bc_index_q=[0 0 0 0];
GDOF_q.P_g=Mesh_node_xy';
GDOF_q.T_g=Mesh_trig_node';
[node_type_q,DOF_q] = global_dof_p_2D_fe_Lagrange(domain,bc_index_q,GDOF_q);

q_gd = FE_u_gd_2D_Lagrange(psi_d_name,node_type_q,GDOF_q);


n_step=1001:8101;
N=513*257;
psi_snap_dave=zeros(N,length(n_step));
u_snap_dave=zeros(N,length(n_step));
v_snap_dave=zeros(N,length(n_step));

S=M_xx+M_yy;
tM=M(DOF_q,DOF_q);
tS=S(DOF_q,DOF_q);
%tM_y0=M_y0(DOF_q,DOF_q);
%tM_x0=M_x0(DOF_q,DOF_q);




dx=1/256;
dy=dx;
xgrid=0:dx:1;
ygrid=0:dy:2;

for k=1:length(n_step)
    
    %--------compute the psi (stream function ) first
    i=n_step(k);
    vor=vor_snap(:,i);
    tvor=vor(DOF_q);
    tpsi=tS\(tM*tvor);
    psi=FE_tu_fe_to_u_fe_2D_Lagrange(tpsi,q_gd,node_type_q,DOF_q);    
    psi_snap_dave(:,k)=psi;
    %psi_map(:,:,i)=reshape(psi,[513 257]);
    %-------Compute the velocity
    %----------------------------
    p=reshape(psi,[513 257]);
    [px, py]=compute_velocity(p,dx,dy,257,513);
    u=py;
    v=-px;
    
    %u_map(:,:,i)=reshape(u, [513 257]);
    %v_map(:,:,i)=reshape(v, [513 257]);
    u_snap_dave(:,k)=reshape(u, [N 1]);
    v_snap_dave(:,k)=reshape(v, [N 1]);
    k,
end