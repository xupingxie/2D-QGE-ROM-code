
function [px, py]=compute_velocity(p,dx,dy,xdim,ydim)

%xdim=length(xgrid);
%ydim=length(ygrid);
% dx=1/256;
% dy=1/256;
% X=0:dx:1;
% Y=0:dy:2;
% 
% idx=1:20:257;
% idy=1:20:513;
% 
% psi=psi_snap(:,1001);
%p=reshape(psi,[513 257]);
px=zeros(ydim,xdim);
py=zeros(ydim,xdim);
% for i =1:ydim
%     for j=1:xdim
%         if i==1
%             py(i,j)=(-(3/2)*p(i,j)+2*p(i+1,j)-(1/2)*p(i+2,j))/dy;
%         elseif i==ydim
%             py(i,j)=((3/2)*p(i,j)-2*p(i-1,j)+(1/2)*p(i-2,j))/dy;
%         else
%             py(i,j)=((1/2)*p(i+1,j)-(1/2)*p(i-1,j))/dy;
%         end
%         
%         if j==1
%             px(i,j)=(-(3/2)*p(i,j)+2*p(i,j+1)-(1/2)*p(i,j+2))/dx;
%         elseif j==xdim
%             px(i,j)=((3/2)*p(i,j)-2*p(i,j-1)+(1/2)*p(i,j-2))/dx;
%         else
%             px(i,j)=((1/2)*p(i,j+1)-(1/2)*p(i,j-1))/dy;
%         end
%     end
% end
%px=zeros(ydim,xdim);

for i =1:ydim
    for j=1:xdim
        if i==1
            py(i,j)=(-25*p(i,j)/12+4*p(i+1,j)-3*p(i+2,j)+4*p(i+3,j)/3-p(i+4,j)/4)/dy;
        elseif i==ydim
            py(i,j)=(25*p(i,j)/12-4*p(i-1,j)+3*p(i-2,j)-4*p(i-3,j)/3+p(i-4,j)/4)/dy;
        elseif i==2
            py(i,j)=(-25*p(i,j)/12+4*p(i+1,j)-3*p(i+2,j)+4*p(i+3,j)/3-p(i+4,j)/4)/dy;
        elseif i==ydim-1
            py(i,j)=(25*p(i,j)/12-4*p(i-1,j)+3*p(i-2,j)-4*p(i-3,j)/3+p(i-4,j)/4)/dy;
        else
            py(i,j)=(p(i-2,j)/12-2*p(i-1,j)/3+2*p(i+1,j)/3-p(i+2,j)/12)/dy;
        end
        
        if j==1
            px(i,j)=(-25*p(i,j)/12+4*p(i,j+1)-3*p(i,j+2)+4*p(i,j+3)/3-p(i,j+4)/4)/dy;
        elseif j==xdim
            px(i,j)=(25*p(i,j)/12-4*p(i,j-1)+3*p(i,j-2)-4*p(i,j-3)/3+p(i,j-4)/4)/dy;
        elseif j==2
            px(i,j)=(-25*p(i,j)/12+4*p(i,j+1)-3*p(i,j+2)+4*p(i,j+3)/3-p(i,j+4)/4)/dy;
        elseif j==xdim-1
            px(i,j)=(25*p(i,j)/12-4*p(i,j-1)+3*p(i,j-2)-4*p(i,j-3)/3+p(i,j-4)/4)/dy;
        else
            px(i,j)=(p(i,j-2)/12-2*p(i,j-1)/3+2*p(i,j+1)/3-p(i,j+2)/12)/dy;
        end
    end
end


% figure(1)
% hold on
% contour(X,Y,p,20);
% quiver(X(idx),Y(idy),py(idy,idx),-px(idy,idx));
% hold off
% axis equal tight
% title('$\psi_{10}$ from FD','interpreter','latex','fontsize',20);
% print('-fillpage','velocity_FD','-dpdf')
% 
% figure(3)
% contourf(X,Y,w_fe_map,20);
% colormap jet;
% axis equal tight;
% title('$\omega_{10}$ from FE projection','interpreter','latex','fontsize',20);
% colorbar;
% print('-fillpage','vorticity_FE','-dpdf')
