function W=simpson_2D_weights(nx,ny,dx,dy)

%xmin;xmax;
%ymin;ymax;

%ny=513;
%nx=257;

W=zeros(ny,nx);


h=dx*dy./9;
for i=1:ny
    for j=1:nx
        
        if i==1 || i==ny || j==1 || j==nx
            
            if i==1 || i==ny
                if mod(j,2)==0
                    W(i,j)=4*h;
                else
                    W(i,j)=2*h;
                end
            end
            
            
            if j==1 || j==nx
                if mod(i,2)==0
                    W(i,j)=4*h;
                else
                    W(i,j)=2*h;
                end
            end
        else
            if mod(i,2)==0 
                if mod(j,2)==0
                    W(i,j)=16*h;
                else
                    W(i,j)=8*h;
                end
            else
                if mod(j,2)==0
                    W(i,j)=8*h;
                else
                    W(i,j)=4*h;
                end
            end
        end
        
    end
end       
                
W(1,1)  = 1*h;   
W(1,nx) = 1*h;
W(ny,1) = 1*h;
W(ny,nx)= 1*h;


% 
% for i=1:nt
%     for j=i:nt
%         
%         p(:,:,i)+p(:,:,j)


        