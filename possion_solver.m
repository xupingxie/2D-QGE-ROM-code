function [v] = possion_solver( hx,Mx,My,f )
%% Mx, number of intervals in y direction 
%% My, # of intervals in x direction
% Mx=512; My=256;

%hx=(Omega(2)-Omega(1))/Mx;
%hy=(Omega(4)-Omega(3))/My;
%x=Omega(1):hx:Omega(2);
%y=Omega(3):hy:Omega(4);

v=zeros(Mx+1,My+1);
MM=(Mx+1)*(My+1);
rhs=zeros(MM,1);
A=zeros(MM,MM);

ne=1/(6*hx^2);
se=ne;
wn=ne;
ws=ne;

cn=2/(3*hx^2);
cs=cn;
ce=cn;
cw=cn;
c=-10/(3*hx^2);

for i=1:Mx+1
    for j=1:My+1
        m=glbidx(i,j,Mx,My);
        if(i==1||i==Mx+1||j==1||j==My+1)
            A(m,m)=1;
            %rhs(m)=feval(u_ext,x(i),y(j));
            rhs(m)=f(i,j);
        else
            n=glbidx(i,j+1,Mx,My);
            e=glbidx(i+1,j,Mx,My);
            s=glbidx(i,j-1,Mx,My);
            w=glbidx(i-1,j,Mx,My);
            
            cne=glbidx(i+1,j+1,Mx,My);
            cse=glbidx(i+1,j-1,Mx,My);
            cwn=glbidx(i-1,j+1,Mx,My);
            cws=glbidx(i-1,j-1,Mx,My);
            
            A(m,m)=c;
            A(m,n)=cn;
            A(m,e)=ce;
            A(m,w)=cw;
            A(m,s)=cs;
            A(m,cne)=ne;
            A(m,cse)=ne;
            A(m,cwn)=ne;
            A(m,cws)=ne;
            
            %rhs(m)=(feval(f_src,x(i-1),y(j))+feval(f_src,x(i+1),y(j))+feval(f_src,x(i),y(j-1))+feval(f_src,x(i),y(j+1))+8*feval(f_src,x(i),y(j)))/12;
            rhs(m)=(f(i-1,j)+f(i+1,j)+f(i,j-1)+f(i,j+1)+8*f(i,j))/12;
        end
    end
    i,
end

vec=A\rhs;

v=reshape(vec,Mx+1,My+1);
    function m=glbidx(i,j,Mx,My)
        m=i+(j-1)*(Mx+1);
    end


end



 

        



