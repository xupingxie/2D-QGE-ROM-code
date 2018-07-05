%%------Date: July 4th 2016

r=10;
Re=450;
Ro=0.0036;


Basis_add_bubble=0;
Num_conforming_nodes_per_element=3;
Load_var_cell=[0,0,1];
t=0;
Load_func=@(y)sin(pi*(y-1));

compute_basis=0;
assemble_tensor=0;
compute_load_vec=0;
compute_snap=0;
directory='Matrix/';
a=0;

%dt=0.01;
time=[10,80];
%t_start=time(1)/dt;
%t_end=time(2)/dt;




if compute_snap==1
    fprintf('load data and sample snapshos\n')
    %load vor_snap_40_50.mat
    skip=10;
    load FE_matrix_top.mat
    t_set=1001:skip:8001;
    w_snap=vor_snap(:,t_set);   %-----sample data for POD
    w_true=vor_snap(:,1001:8001);
    save(strcat(directory,'snapdata_','time_',num2str(time(1)),'_',num2str(time(2))),'w_snap','w_true','-v7.3')
     
 else
     fprintf('Load data FEmesh and snap\n')
     load FE_matrix_bot.mat
     %load(strcat(directory,'snapdata_','time_',num2str(time(1)),'_',num2str(time(2))))
 end



elem=Mesh_trig_node;
nodes=Mesh_node_xy;
FEM.elem=elem;
FEM.nodes=nodes;
S=M_xx+M_yy;

%%%---------------compute the load vector
if compute_load_vec==1
    fprintf('compute the source load\n')
    FE_load=FEM_load_2D(Load_func, Load_var_cell, nodes, elem, t,...
        Basis_add_bubble, Num_conforming_nodes_per_element);
    save(strcat(directory,'FE_load_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2))),'FE_load')
else
    fprintf('load vector\n')
    load(strcat(directory,'FE_load_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2))),'FE_load')
end


%%%------------compute the POD basis

if compute_basis==1
    if a==0
        fprintf('Compute the POD basis\n')
        u_snap=w_snap;
        [POD_vor,POD_psi]=compute_POD_w_psi_basis(u_snap,M,S,r,a,FEM);
    else
        fprintf('Compute the FTLE-POD basis\n')
        load ftle_snap.mat
        u_snap=[w_snap;sqrt(a)*ftle_snap];
        [POD_vor,POD_psi]=compute_POD_w_psi_basis(u_snap,M,S,r,a,FEM);
    end
    save(strcat(directory,'POD_data_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2)),'_a',num2str(a)),'POD_vor','POD_psi')
else
    fprintf('Load POD basis\n')
    load(strcat(directory,'POD_data_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2)),'_a',num2str(a)),'POD_vor','POD_psi')
end
fprintf('Complete POD basis!\n')


%%%%------------------Compute Tesnor
if assemble_tensor==1
    fprintf('Assemble the Tensor\n')
    [Txy, Tyx]=POD_tensor_assemble2D_DG(POD_vor,nodes, elem);
    save(strcat(directory,'tensor_data_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2)),'_a',num2str(a)),'Txy','Tyx');
    %[T]=POD_tensor_assemble2D_new_QGE2(POD_vor,nodes, elem);
    %save(strcat(directory,'tensor_data_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2))),'T');
else
    fprintf('load tensor\n')
    %[Tx, Ty]=POD_tensor_assemble2D_DG(POD_vor,nodes, elem)
    %load(strcat(directory,'tensor2_data_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2))));
    load(strcat(directory,'tensor_data_r',num2str(r),'time_',num2str(time(1)),'_',num2str(time(2)),'_a',num2str(a)))
end


%[B,P,Q]=QGE_ROM_matrices_assemble_nomean2(POD_vor,POD_psi,Txy,Tyx,M_x0,S,FE_load,Re,Ro);
%[B,P,Q]=QGE_ROM_matrices_assemble_nomean(POD_vor,POD_psi,T,M_x0,S,FE_load,Re,Ro);
w_POD=POD_vor;
p_POD=POD_psi;
Sr = (1/Re)*w_POD'*S*w_POD;
Mr_x0 = (1/Ro)*p_POD'*M_x0*w_POD;
fr = (1/Ro)*w_POD'*FE_load;
lmM = POD_vor'*M*
Tenx = zeros(r,r,r);
Teny = zeros(r,r,r);
%for i=1:r
%    Tenx(:,:,i) = p_POD'*Txy(i).ux*w_POD;        %(\varphi_psi_x,varphi_p_y,\varphi)
%    Teny(:,:,i) = p_POD'*Tyx(i).uy*w_POD;        %(\varphi_psi_y,varphi_p_x,\varphi)
%    %Tenx(:,:,i)=p_POD'*T(i).xy*w_POD;
%    %Teny(:,:,i)=p_POD'*T(i).yx*w_POD;
%end
for i=1:r
    Tenx(:,:,i) = w_POD'*Txy(i).ux*p_POD;        %(\varphi_psi_x,varphi_p_y,\varphi)
    Teny(:,:,i) = w_POD'*Tyx(i).uy*p_POD;        %(\varphi_psi_y,varphi_p_x,\varphi)
    %Tenx(:,:,i)=p_POD'*T(i).xy*w_POD;
    %Teny(:,:,i)=p_POD'*T(i).yx*w_POD;
end
fprintf('Complete offline matrix assembling \n ')

fprintf('computing ODE now\n')
ddt=0.001;
n_time=(time(2)-time(1))/ddt;
save_index = 10;
t0=time(1);
wr0=POD_vor'*M*w_snap(:,1);
Cw = zeros(r,7001);
Cw(:,1)=wr0;

for k=1:n_time
    
    %--------Compute K1
    t=t0;
    a=wr0;
    %K1=Assemble_differential_F(a,r,fr1,P,Q);
    NonL=nonlinear_term(a,Tenx,Teny,r);
    K1=-Sr*a+Mr_x0*a+NonL+fr;
    
    %--------Compute K2
    t=t0+ddt/2;
    a=wr0+ddt*K1/2;
    %K2=Assemble_differential_F(a,r,fr2,P,Q);
    NonL=nonlinear_term(a,Tenx,Teny,r);
    K2=-Sr*a+Mr_x0*a+NonL+fr;
    
    %---------Compute K3
    a=wr0+ddt*K2/2;
    %K3=Assemble_differential_F(a,r,fr2,P,Q);
    NonL=nonlinear_term(a,Tenx,Teny,r);
    K3=-Sr*a+Mr_x0*a+NonL+fr;
    
    %---------Compute K4
    t=t0+ddt;
    a=wr0+ddt*K3;
    %K4=Assemble_differential_F(a,r,fr3,P,Q);
    NonL=nonlinear_term(a,Tenx,Teny,r);
    K4=-Sr*a+Mr_x0*a+NonL+fr;
    
    %----------Compute w^n+1
    
    w_new=wr0+ddt*(K1+2*K2+2*K3+K4)/6;
    
    if rem(k,save_index) == 0
        Cw(:,k/save_index+1) = w_new;
        k,t,
        %Cw(:,k+1)=C1;
    end
    wr0=w_new;
    t0=t;
end
%options_ODE=odeset('Reltol', 1e-8,'AbsTol', 1e-8);
%y0=POD_vor'*M*w_snap(:,1);
%[~,temp]=ode45(@(t,y)Galerkin_QGE(y,r,B,P,Q),(time(1):0.01:time(2)),y0,options_ODE);

w_rom=POD_vor*Cw;
psi_rom=POD_psi*Cw;


%%%%%%------compute velocity
% fprintf('Compute the velocity field\n')
% n_step=size(psi_rom,2);
% for i=1:n_step
%     
%     %--------compute the vorticity ( function ) first
%     psi=psi_rom(:,i);
%     %tpsi=psi(DOF_q);
%     %tvor=tM\(tS*tpsi);
%     %vor=FE_tu_fe_to_u_fe_2D_Lagrange(tvor,q_gd,node_type_q,DOF_q);    
%     %vor_snap(:,i)=vor;
%     %psi_map(:,:,i)=reshape(psi,[513 257]);
%     %-------Compute the velocity
%     %----------------------------
% 
%     u=M\(M_y0*psi);
%     v=M\(-M_x0*psi);
% %     
% %     %u_map(:,:,i)=reshape(u, [513 257]);
% %     %v_map(:,:,i)=reshape(v, [513 257]);
%      u_grom(:,i)=u;
%      v_grom(:,i)=v;
%     i,
% end
% 
% fprintf('Done!\n')



% for i=1:1000
%     h=pcolor(F(:,:,i));set(h,'edgecolor','none');colormap jet;
%     pause(.1)
% end
% 
% N=513*257;
% ftle_snap=zeros(131841,1001);
% for i=1:1000
%     ftle_snap(:,i)=reshape(F(:,:,i),[N 1]);
% end
% u_dns=zeros(131841,1001);
% v_dns=zeros(131841,1001);
% for i=1:n_step
%     
%     %--------compute the vorticity ( function ) first
%     psi=psi_snap(:,i);
% 
%     u=M\(M_y0*psi);
%     v=M\(-M_x0*psi);
% 
%     u_dns(:,i)=u;
%     v_dns(:,i)=v;
%     i,
% end
E=[];
for i=1:7001
    err=(psi_snap_10_81(:,i)-psi_rom(:,i))'*M*(psi_snap_10_81(:,i)-psi_rom(:,i));
    E=[E,err];
end
EL2=sum(E)./7001;
% E=[];
% for i=1:7001;
%     er=w_true(:,i)'*M*w_true(:,i);
%     E=[E,er];
% end
% Er=sum(E)./7001;
% 
% E=[];
% for i=1:7001;
%     ftle=reshape(F(:,:,i),[131841 1]);
%     er=ftle'*M*ftle;
%     E=[E,er];
% end
% Er=sum(E)./7001;
% 
% for i=1:7001;
%     ww=ww+F(:,:,i);
% end
% wmean=ww./7001;