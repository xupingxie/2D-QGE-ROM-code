function dy=Galerkin_QGE(y,r,B,P,Txy,Tyx)

dy=zeros(r,1);
%-psi_POD'*T(i).yx*w_POD+psi_POD'*T(i).xy*w_POD;
for i=1:r
    dy(i)=B(i,1)+P(i,:)*y-y'*Tyx(:,:,i)*y+y'*Txy(:,:,i)*y;
end

end