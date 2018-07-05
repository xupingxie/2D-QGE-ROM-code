function [node_type, DOF_p]=global_dof_p_2D_fe_Lagrange(domain,bc_index,GDOF)

tol=eps;
node_type=zeros(2,size(GDOF.P_g,2));
%node_type=-200*ones(2,size(GDOF.P_g,2));
DOF_p=zeros(1,size(GDOF.P_g,2));
unknown_count=0;

for i=1:length(GDOF.P_g)
    if abs(GDOF.P_g(2,i)-domain(3))<=tol
        if (GDOF.P_g(1,i)==domain(1)&&GDOF.P_g(2,i)==domain(3))
            if bc_index(1)~=bc_index(4)
                node_type(2,i)=0;
            else
                node_type(2,i)=bc_index(1);
            end
        elseif GDOF.P_g(1,i)==domain(2)&&GDOF.P_g(2,i)==domain(3)
            if bc_index(1)~=bc_index(4)
                node_type(2,i)=0;
            else
                node_type(2,i)=bc_index(1);
            end
        else       
            node_type(2,i)=bc_index(1);
        end
        if node_type(2,i)<0
              unknown_count=unknown_count+1;
              node_type(1,i)=unknown_count;
              DOF_p(unknown_count)=i;
        end
    elseif abs(GDOF.P_g(1,i)-domain(2))<=tol     
        node_type(2,i)=bc_index(2);
        if bc_index(2)<0
            unknown_count=unknown_count+1;
            node_type(1,i)=unknown_count;
            DOF_p(unknown_count)=i;
        end
    elseif abs(GDOF.P_g(2,i)-domain(4))<=tol
        if (GDOF.P_g(1,i)==domain(1)&&GDOF.P_g(2,i)==domain(4))
            if bc_index(3)~=bc_index(4)
                node_type(2,i)=0;
            else
                node_type(2,i)=bc_index(3);
            end
        elseif GDOF.P_g(1,i)==domain(2)&&GDOF.P_g(2,i)==domain(4)
            if bc_index(1)~=bc_index(4)
                node_type(2,i)=0;
            else
                node_type(2,i)=bc_index(3);
            end
        else       
            node_type(2,i)=bc_index(3);
        end
        if bc_index(3)<0
            unknown_count=unknown_count+1;
            node_type(1,i)=unknown_count;
            DOF_p(unknown_count)=i;
        end
    elseif abs(GDOF.P_g(1,i)-domain(1))<=tol
        node_type(2,i)=bc_index(4);
        if bc_index(4)<0
            unknown_count=unknown_count+1;
            node_type(1,i)=unknown_count;
            DOF_p(unknown_count)=i;
        end
    else
         node_type(2,i)=-200;
         unknown_count=unknown_count+1;
         node_type(1,i)=unknown_count;
         DOF_p(unknown_count)=i;
    end
end

DOF_p=DOF_p(1:unknown_count);
return;

end