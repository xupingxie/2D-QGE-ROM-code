function [POD,Diag_S,d,CumEng,CumEng_ratio,rows,dim] = ...
    PODbasis_QGE(u_snap,M)
%============================================================================
% generate POD basis by the method of snapshots
% choose different eigen problems basing on the shape of u: tall or low
% the weight matrix is M
%============================================================================
if nargin<6
    fig_num = 5;
end
[nm,~] = size(M);
[n ,m] = size(u_snap);
dim    = n/nm;
% if dim is the dimension of velocity, in incompressible flow,
% do pod on [u,v,w] together, to keep the divergence free property.
rows   = 1:nm:n+1;

u_snap      = 1/sqrt(m)*u_snap;   %----Zhu's approach

fprintf(1,'Solve POD by method of snapshots\n');
for i = 1:dim
    uc = u_snap(rows(i):rows(i+1)-1,:);
    if i ==1
        Energy = uc'*M*uc;
    else
        Energy = Energy + uc'*M*uc;
    end
end




d      = rank(full(Energy));
fprintf(1,['Rank of data set is ',num2str(d),'\n']);
POD    = zeros(n, d);
%u_POD = zeros(n/dim,d);

%ftle_POD =u_POD;

if d<size(Energy,1)-1
    [V,Lambda] = eigs(Energy,d);
else
    [V,Lambda] = eig(Energy);
end
Diag_S = sqrt(diag(Lambda(1:d, 1:d)));
%Diag_S = sqrt(diag(Lambda));
index=1;
while index~=0
    Num_change_index=nnz(Diag_S - sort(Diag_S,'descend')); % check if eigenvalues in descending order
    if Num_change_index==0
        index = 0;
    else
        fprintf(1,[num2str(index),' not in descending order\n']);
        index=index+1;
        [Diag_S, Ind]=sort(Diag_S,'descend');
        V = V(:, Ind);
        change_index=find(Diag_S - sort(Diag_S,'descend')~=0);
        if isempty(change_index)==0
            for i=1:2:Num_change_index % in this special case, wrong order always appears in couples
                temp_1 = change_index(i);
                temp_2 = change_index(i+1);
                temp_3 = Diag_S(temp_2);         temp4       = V(:,temp_2);
                Diag_S(temp_2) = Diag_S(temp_1); V(:,temp_2) = V(:,temp_1);
                Diag_S(temp_1) = temp_3;         V(:,temp_1) =temp4;
            end
        end
    end
end
%     for i=1:dim
%         for k=1:d
%             POD(rows(i):rows(i+1)-1,k) = 1/Diag_S(k)*u(rows(i):rows(i+1)-1, :)*V(:,k);
%         end
%     end
for i=1:dim
    temp = u_snap(rows(i):rows(i+1)-1, :)*V;
    POD(rows(i):rows(i+1)-1,:) = temp./(ones(nm,1)*Diag_S(1:d,1)');
    %POD(rows(i):rows(i+1)-1,:) = temp./(ones(nm,1)*Diag_S');
end

%u_POD = POD(rows(1):rows(2)-1,:);
%v_POD = POD(rows(2):rows(3)-1,:);
%ftle_POD = POD(rows(3):rows(4)-1,:);

% fprintf(1, 'Note Diag_S is sqrt{lambda}\n');
 CumEng = zeros(1,d);
 for i=1:d
     CumEng(i) = sum(Diag_S(1:i).^2);
 end
 CumEng_ratio = CumEng/sum(Diag_S.^2);
fprintf(1, 'The first 10 POD basis functions give the energy share as follows.\n');
disp([(1:10)', CumEng_ratio(1:10)']);



%figure
%subplot(1,2,1), semilogy(1:d, Diag_S(1:d).^2,'r*');
%subplot(1,2,2), plot(1:d,CumEng_ratio,'r*')
%Fig_name = [fig_pod,'_eig_eng_distribution'];
%print('-f', '-depsc',Fig_name)
%fig_num  = fig_num + 1;
%semilogy(1:d,Diag_S(1:d).^2,'o-');
%xlabel('r')
%ylabel('\lambda_i / \Sigma \lambda_i')
%grid on
%title('Eigenvalues')


end


