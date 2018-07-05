function [T]=POD_tensor_assemble2D_new_QGE(POD_vor,x, e_conn)
%---- this is the new version 
%---- June, 7, 2016

% To evaluate the trilinear term

% ouput T_yx=(h_y phi_k,h_x); T_xy=(h_x phi_k,h_x);

[N, r]=size(POD_vor);

%[N,          dim    ] = size(x     );
[n_elements, nel_dof] = size(e_conn);

%Tux(q).loc = []; Tx(q).ux = [];XT_x(q).ux=[];
%Tuy(q).loc = []; Ty(q).uy = [];XT_y(q).uy=[];

% Due to memory limitations, we will fill the mass matrix by integrating
% over partitions of elements.  The constant (max_elem_per_partition)
% should be adjusted based on available memory.
max_elem_per_partition = 500;

n_part = floor( n_elements / (max_elem_per_partition+1) ) + 1;
elem_segment = floor( linspace(0,n_elements,n_part+1) );
max_part_size = max( diff( elem_segment ) );


    [rr,ss,ww] = FEM_twod_gauss(7);
    %one = ones(size(w));
    
    for n_pt = 1:n_part
        II   = sparse( max_part_size*nel_dof^2,1 );
        JJ   = II;
   
        for i=1:r
            XT(i).xy = sparse(max_part_size*nel_dof^2,1);
            XT(i).yx = sparse(max_part_size*nel_dof^2,1);
        end
        entry_counter = 0;
        for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
            nodes_local             = e_conn(n_el,:);
            x_local                 = x(nodes_local,:);
            [~, w_g, phi, p_x, p_y] = FEM_twod_shape(x_local,rr,ss,ww);
            
            %   Assemble into global matrix
            %-----------------
            
            for i=1:nel_dof
                for j=1:nel_dof
                    entry_counter=entry_counter+1;
                    II(entry_counter) = nodes_local(i);
                    JJ(entry_counter) = nodes_local(j);
                    
                    for k=1:r
                        podu_k = phi*POD_vor(nodes_local,k);
                        
                        Txy_local=FEM_twod_bilinear(podu_k, p_x, p_y, w_g); 
                        Tyx_local=FEM_twod_bilinear(podu_k, p_y, p_x, w_g); 
                        
                        XT(k).xy(entry_counter) = iszero(Txy_local(i,j));
                        XT(k).yx(entry_counter) = iszero(Tyx_local(i,j));
                        %Tux(k).loc = FEM_twod_bilinear(podu_k, p_x, p_y, w_g);    %(\varphi_x()*h_y, h)
                        %Tuy(k).loc = FEM_twod_bilinear(podu_k, p_y, p_x, w_g);    %(\varphi_y()*h_x, h)
                    end
                end
            end
            
            if mod(n_el,50) == 0 && mod(n_el, 500)~=0
                fprintf(1, [num2str(n_el),'\t']);
            elseif mod(n_el,500) ==0
                fprintf(1, [num2str(n_el),'\n']);
            end
        end
            
        if ( n_pt==1 )
            
            for k=1:r
                T(k).xy = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XT(k).xy(1:entry_counter),...
                    N, N );
                T(k).yx = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XT(k).yx(1:entry_counter),...
                    N, N );
            end
            
        else
            for k=1:r
                T(k).xy = T(k).xy + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XT(k).xy(1:entry_counter),...
                    N, N );
                T(k).yx = T(k).yx + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XT(k).yx(1:entry_counter),...
                    N, N );
            end
        end
        
    end
            
    function entry=iszero(entry)
        if abs(entry) < eps
            entry=0;
        end