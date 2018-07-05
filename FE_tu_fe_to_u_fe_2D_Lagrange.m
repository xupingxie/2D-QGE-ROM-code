function u_fe=FE_tu_fe_to_u_fe_2D_Lagrange(tu_fe,u_gd,node_type,DOF_p)

N=size(node_type,2);
u_fe=zeros(N,1);
tol=eps;

u_fe(DOF_p,1)=tu_fe;

for i=1:N
    if abs(node_type(2,i)-0)<tol
        u_fe(i)=u_gd(i);
    end
end

        
        








