function [F,P,Q]=QGE_ROM_matrices_assemble_nomean(POD_vor,POD_psi,T,M_x0,S,FE_load,Re,Ro)

%%---Assemble Matrices,
[N, r]= size(POD_vor);

%psi_bar=psi_mean;
%w_bar=vor_mean;
w_POD=POD_vor;
psi_POD=POD_psi;

%for i=1:r
    %L(i,1)=-w_bar'*S*w_POD(:,i);
    %N(i,1)=-psi_bar'*T(i).yx*w_bar+psi_bar'*T(i).xy*w_bar;
    %H(i,1)=psi_bar'*M_x0*w_POD(:,i);
    %F(i,1)=w_POD(:,i)'*FE_load;
    %B(i,1)=(1/Re)*L(i,1)+N(i,1)+(1/Ro)*H(i,1)+F(i,1);
%    F(i,1)=(1/Ro)*FE_load'*w_POD(:,i);
%end
F=(1/Ro)*w_POD'*FE_load;
%B=F;

P=zeros(r,r);
%P_temp=zeros(r,r);
for i=1:r
    %mat_L(i,:)=-w_POD'*S*w_POD;
    %mat_N1(i,:)=-psi_bar'*T(i).yx*w_POD+psi_bar'*T(i).xy*w_POD;
    %mat_N2(:,i)=-psi_POD'*T(i).yx*w_bar+psi_POD'*T(i).xy*w_bar;
    %mat_N(i,:)=mat_N1(i,:)+mat_N2(:,i)';
    %mat_H(i,:)=psi_POD'*M_x0*w_POD(:,i);
    %P(i,:)=(1/Re)*mat_L(i,:)+mat_N(i,:)+(1/Ro)*mat_H(i,:);
    %P_temp(i,:)=mat_N(i,:)+(1/Ro)*mat_H(i,:);
    
    P1(i,:)=-(1/Re)*w_POD'*S*w_POD(:,i);
    P2(i,:)=(1/Ro)*psi_POD'*M_x0*w_POD(:,i);
    P(i,:)=P1(i,:)+P2(i,:);
    %P_temp(i,:)=(1/Ro)*psi_POD'*M_x0*w_POD(:,i);
end
%P=P_temp-(1/Re)*w_POD'*S*w_POD;

Q=cell(r,1);
for i=1:r
    %T(i).yx=Tyx(i).uy;
    %T(i).xy=Txy(i).ux;
    %T(i).yx=Tyx(i).uy;
    %T(i).xy=Txy(i).ux;
    
    Q{i}(:,:)=-psi_POD'*T(i).yx*w_POD+psi_POD'*T(i).xy*w_POD;
    %Q{i}(:,:)=-psi_POD'*T(i).xy*w_POD+psi_POD'*T(i).yx*w_POD;
end
