function N=nonlinear_term(C0,Tenx,Teny,r)

    Nh1 = zeros(r,1);
    Nh2 = zeros(r,1);
    for i=1:r
        Nh1(i) = C0'*Tenx(:,:,i)*C0;
        Nh2(i) = C0'*Teny(:,:,i)*C0;
    end
    
    N=-Nh2+Nh1;
    %N=-Nh1+Nh2;
end

