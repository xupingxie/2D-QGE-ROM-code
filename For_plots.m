N=513*257;
ftle_snap=zeros(N,7001);

for i=1:7001
    w=reshape(F(:,:,i),[N 1]);
    ftle_snap(:,i)=w;
end


figure(1)
contourf(xgrid,ygrid,reshape(psimean,[513 257]),20);colormap jet;axis equal tight;
title('mean $\psi$, $r=30$','interpreter','latex','fontsize',20)
color bar
figure(1)
k=2001;
contourf(xgrid,ygrid,reshape(psi_rom(:,k),[513 257]),20);colormap jet;axis equal tight;
title('$\psi(t=60)$, $r=30$','interpreter','latex','fontsize',20)



figure(1)
k=6501; %k=2001;5001
h=pcolor(xgrid,ygrid,F(:,:,k));
set(h,'edgecolor','none');
colormap jet;
axis equal tight
title('FTLE $\lambda$ at $t=75$','interpreter','latex','fontsize',20)

figure(1)
k=5001; %k=2001;5001
h=pcolor(xgrid,ygrid,F(:,:,k));
set(h,'edgecolor','none');
colormap jet;
axis equal tight
title('$\lambda(t=60)$ ','interpreter','latex','fontsize',20)

figure(2)
imagesc(xgrid,ygrid,F(:,:,k));
colormap jet;
axis equal tight

figure(3)
a=F(:,:,k);
xx=3:510;yy=2:254;
h=surf(X(xx,yy),Y(xx,yy),a(xx,yy));
set(h,'LineStyle','none')
view([15 -20 38])
%colormap jet;

figure(4)
mesh(X,Y,F(:,:,k));

t=10:0.01:80;
hold on
loglog(t,Energy_dns,'b');
loglog(t,Energy_ROMr10,'r');
loglog(t,Energy_ROMr10_f,'k');
loglog(t,Energy_ROMr10_a1,'c');
hold off
hlegend=legend('DNS','GROM','$\alpha-ROM$','$\lambda-ROM$');
set(hlegend,'interpreter','latex','Fontsize',17)
title('Time evolution of KE','fontsize',17);
xlabel('t','fontsize',17)
ylabel('E(t)','fontsize',17)

skip=10;
idx=1:skip:257;
idy=1:skip:513;
p=reshape(psi_rom(:,end),[513 257]);
u1=u_snap_meanftle3(:,:,end);
v1=v_snap_meanftle3(:,:,end);

hold onn
contour(xgrid,ygrid,p,20);
quiver(xgrid(idx),ygrid(idy),u1(idy,idx),v1(idy,idx));
hold off

N=513*257;
u_meanftle1=zeros(N,2001);
u_meanftle2=zeros(N,2000);
u_meanftle3=zeros(N,3100);
for i=1:2001
    u_meanftle1(:,i)=reshape(u_snap_meanftle1(:,:,i),[N 1]);
end

for i=2002:4001
    u_meanftle2(:,i)=reshape(u_snap_meanftle2(:,:,i),[N 1]);
end
for i=1:3001
    u_meanftle3(:,i)=reshape(u_snap_meanftle3(:,:,i),[N 1]);
end

v_meanftle1=zeros(N,2001);
v_meanftle2=zeros(N,2000);
v_meanftle3=zeros(N,3100);
for i=1:2001
    v_meanftle1(:,i)=reshape(v_snap_meanftle1(:,:,i),[N 1]);
end

for i=1:2000
    v_meanftle2(:,i)=reshape(v_snap_meanftle2(:,:,i),[N 1]);
end
for i=1:3001
    v_meanftle3(:,i)=reshape(v_snap_meanftle3(:,:,i),[N 1]);
end

N=513*257;
a1ftle_r30=zeros(N,7001);
for i=1:7001
    a1ftle_r30(:,i)=reshape(F(:,:,i),[N 1]);
end

pmean=mean(psi_rom,2);
figure(1)
xgrid=0:1/256:1;
ygrid=0:1/256:2;
%p=reshape(POD_psi(:,5),[513 257]);
%p=F(:,:,3001);
p=reshape(pmean,[513 257]);
contourf(xgrid,ygrid,p,20);
colormap jet;
axis equal tight

pmean=POD_psi(:,10);
figure(1)
xgrid=0:1/256:1;
ygrid=0:1/256:2;
p=reshape(pmean,[513 257]);
contourf(xgrid,ygrid,p,20);
colormap jet;
axis equal tight

t=10:0.01:80;
hold on
plot(t,Energy_dns,'b');
plot(t,Energy_ROMr10,'r');
plot(t,Energy_ROMr10_f,'k');
plot(t,Energy_ROMr10_a1,'c');
hold off
hlegend=legend('DNS','GROM','$\alpha-ROM$','$\lambda-ROM$');
set(hlegend,'interpreter','latex','Fontsize',17,'linewidth',2)
title('Time evolution of KE','fontsize',17);
xlabel('t','fontsize',17)
ylabel('E(t)','fontsize',17)
grid on

t=10:0.01:80;
hold on
plot(t,Error_grom_r10,'b');
plot(t,Error_a1rom_r10,'r');
plot(t,Error_from_r10,'k');
hold off
hlegend=legend('GROM','$\alpha-ROM$','$\lambda-ROM$');
set(hlegend,'interpreter','latex','Fontsize',17)
title('Time evolution of error','fontsize',17);
xlabel('t','fontsize',17)
ylabel('$\|e\|_{L^2}$','interpreter','latex','fontsize',17)
grid on

t=10:0.01:80;
hold on
loglog(t,Error_grom_r10,'b');
loglog(t,Error_a1rom_r10,'r');
loglog(t,Error_from_r10,'k');
hold off
hlegend=legend('GROM','$\alpha-ROM$','$\lambda-ROM$');
set(hlegend,'interpreter','latex','Fontsize',17)
title('Time evolution of error','fontsize',17);
xlabel('t','fontsize',17)
ylabel('$\|e\|_{L^2}$','interpreter','latex','fontsize',17)
grid on

N=513*257;
ftle_snap=zeros(N,7001);
for i=1:7001
    p=reshape(F(:,:,i),[N 1]);
    ftle_snap(:,i)=p;
    i,
end

loglog(D,'-*');
title('Eigenvalues of the correlation matrix $\alpha=10^4$','interpreter','latex','fontsize',17)
grid on
%xlabel('$r$','interpreter','latex', 'fontsize',17)

figure(1)
p1=psi_rom(:,1001);
p2=psi_snap_fe(:,1001);
contourf(xgrid,ygrid,reshape(p1,[513 257]),20);
colormap jet;
title('ROM')
figure(2)
contourf(xgrid,ygrid,reshape(p2,[513 257]),20);
colormap jet;
title('DNS')

for i=1:7001
    a=reshape(ftle_true(:,i),[513 257]);
    figure
    h=pcolor(xgrid,ygrid,a);
    set(h,'edgecolor','none');
    colormap jet;
    axis equal tight;
    pause(.1);
end

R=[10, 15, 20, 25, 30];
E_from=[];
for i=1:length(R)
    r=R(i);
    ener=sum(D(1:r))/sum(D);
    E_from=[E_from,ener];
end
Ener=zeros(9,5);

Ener(1,:)=E_arom0;
Ener(2,:)=E_arom1;
Ener(3,:)=E_arom10;
Ener(4,:)=E_arom100;
Ener(5,:)=E_arom1000;
Ener(6,:)=E_arom10000;
Ener(7,:)=E_arom105;
Ener(8,:)=E_arom106;
Ener(9,:)=E_from;

a=[0,1,10,100,10^3,10^4,10^5,10^6];
x=logspace(a);
figure(1)
semilogx(a, Ener(1:8,1),'-*');
hold all
semilogx(a, Ener(1:8,2),'-d');
semilogx(a, Ener(1:8,3),'-s');
semilogx(a, Ener(1:8,4),'-o');
semilogx(a, Ener(1:8,5),'-.');
hlegend=legend('r=10','r=15','r=20','r=25','r=30');
set(hlegend,'Fontsize',14);
title('Energy ratio growing for $\alpha$ inner product','interpreter','latex','fontsize',14);
xlabel('$\alpha$','interpreter','latex','fontsize',14);
ylabel('$\frac{\sum_{j=1}^r\lambda_j}{\sum_{j=1}^d\lambda_j}$','interpreter','latex','fontsize',20);


Error_mean=zeros(5,6);
Error_L2=zeros(5,6);

Error_L2(1,:)= [4.2407e+03, 4.7738e+03, 9.7892e+03, 3.9602e+03, 33.211, 40.4220];
Error_L2(2,:)= [1.8787e+03, 1.6912e+03, 55.5664, 60.5894, 36.7646, 29.5179];
Error_L2(3,:)= [47.5402, 55.1052, 42.6515, 34.6073, 34.7291, 32.1291];
Error_L2(4,:)= [38.8147, 37.2375, 34.9139, 34.2630, 29.8403, 29.4772];
Error_L2(5,:)= [33.7272, 35.0747, 35.6443, 28.8374, 27.2999, 27.6152];

Error_mean(1,:)= [332.9672, 938.7975, 6.5470e+03, 499.0825, 2.4506, 0.8803];
Error_mean(2,:)= [1.0406e+03, 842.7060, 11.1189, 10.6121, 1.1724, 1.9613];
Error_mean(3,:)= [10.0911, 12.1001, 6.8204, 2.5756, 1.1632, 1.2657];
Error_mean(4,:)= [3.0556, 3.5823, 3.0208, 1.8140, 0.4542, 0.3266];
Error_mean(5,:)= [2.8330, 2.1385, 2.2526, 1.1904, 0.4861, 0.5860];
r=[10,15,20,25,30];

figure(1)
semilogy(r,Error_mean(:,1),'--');
hold all
semilogy(r,Error_mean(:,2),'-s');
semilogy(r,Error_mean(:,3),'-*');
semilogy(r,Error_mean(:,4),'-o');
semilogy(r,Error_mean(:,5),'-d');
semilogy(r,Error_mean(:,6),'-.');
hlegend=legend('$\alpha=0$','$\alpha=1$','$\alpha=10$','$\alpha=10^2$','$\alpha=10^3$','$\alpha=10^4$');
set(hlegend,'interpreter','latex','Fontsize',14);
xlabel('$r$','interpreter','latex','fontsize',14);
ylabel('$\|\overline{\psi}-\overline{\psi_r}\|$','interpreter','latex','fontsize',20);
title('The error of mean stream function','fontsize',14);

figure(1)
semilogy(r,Error_L2(:,1),'--');
hold all
semilogy(r,Error_L2(:,2),'-s');
semilogy(r,Error_L2(:,3),'-*');
semilogy(r,Error_L2(:,4),'-o');
semilogy(r,Error_L2(:,5),'-d');
semilogy(r,Error_L2(:,6),'-.');
hlegend=legend('$\alpha=0$','$\alpha=1$','$\alpha=10$','$\alpha=10^2$','$\alpha=10^3$','$\alpha=10^4$');
set(hlegend,'interpreter','latex','Fontsize',14);
xlabel('$r$','interpreter','latex','fontsize',14);
ylabel('$\frac{1}{M}\sum_{j=1}^M\|e_j\|$','interpreter','latex','fontsize',20);
title('The time averaged L2 error','fontsize',14);

t=10:0.01:80;
figure(1)
semilogy(t,Error_arom0_r30);
hold all
semilogy(t,Error_arom1_r30);
semilogy(t,Error_arom10_r30);
semilogy(t,Error_arom100_r30);
semilogy(t,Error_arom1000_r30);
semilogy(t,Error_arom104_r30);
%semilogy(r,Error_L2(:,6),'-.');
hlegend=legend('$\alpha=0$','$\alpha=1$','$\alpha=10$','$\alpha=10^2$','$\alpha=10^3$','$\alpha=10^4$');
set(hlegend,'interpreter','latex','Fontsize',14);
xlabel('$t$','interpreter','latex','fontsize',14);
ylabel('$\|e_j\|_{L^2}$','interpreter','latex','fontsize',20);
title('The time evolution of error r=30','fontsize',14);

figure(1)
semilogy(t,Energy_aROM0_r30);
hold all
semilogy(t,Energy_aROM1_r30);
semilogy(t,Energy_aROM10_r30);
semilogy(t,Energy_aROM100_r30);
semilogy(t,Energy_aROM1000_r30);
semilogy(t,Energy_aROM104_r30);
semilogy(t,Energy_dns,'b');
%semilogy(r,Error_L2(:,6),'-.');
hlegend=legend('$\alpha=0$','$\alpha=1$','$\alpha=10$','$\alpha=10^2$','$\alpha=10^3$','$\alpha=10^4$','DNS');
set(hlegend,'interpreter','latex','Fontsize',12);
xlabel('$t$','interpreter','latex','fontsize',14);
ylabel('$E(t)$','interpreter','latex','fontsize',20);
title('The time evolution of error r=20','fontsize',14);


A=[94,89,92,85,90,60,77,77,71,93,83,74,91,81,86,90,94,81,87,91,53,67,65,70,76,83,77,79,87,90,89];

A=[91,82,74,71,82,50,62,69,66,89,70,62,90,68,77,85,84,70,75,85,46,64,63,62,67,66,72,74,73,82,84];

a=F(:,:,4000);
figure
h=pcolor(xgrid,ygrid,a);
set(h,'edgecolor','none');
colormap jet;
axis equal tight;

a=reshape(ftle_mean, [513 257]);
figure
contourf(xgrid,ygrid,a,20);
colormap jet;
axis equal tight;

chox=10:247;
choy=10:503;
xx=xgrid(chox);
yy=ygrid(choy);
ap=p(choy,chox);
h=pcolor(xx,yy,ap);
set(h,'edgecolor','none');
colormap jet;
axis equal tight;




    



    
    
    
    
    
    
