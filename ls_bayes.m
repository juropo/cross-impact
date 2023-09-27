function [x, C0, resnorm,residual,exitflag,output,lambda]=ls_bayes(C_I,m,scenp,states)
resnorm=0;
residual=0;
if isempty(C_I) || isempty(m) || isempty(scenp) || length(states)<2
    x=m(:);
    C0=[];
    return
end

tic

C0=sqrt(2.^C_I);

N=states(1);
nparents=length(states)-1;


%foo=m*2.^C;

%rajoitusmatriisi
constraints=[];
equals=[];

skenprob2=scenp;

%marginaalirajoitteet

for mi=1:N
    temp=zeros(1,N);
    temp(mi)=1;
    temp2=(skenprob2(:)*temp)';
    constraints=[constraints;temp2(:)'];
    equals=[equals;m(mi)];
end

% Tämä rajoite pitää huolen, että kaikki ehdolliset jakauvat summautuvat
% yhteen, eli sum(cdist(:,i,j,k)=1
% for p1=1:3
%     for p2=1:3
%         for p3=1:3
%             temp2=indexmagic2(2,p1,states).*indexmagic2(3,p2,states)...
%                 .*indexmagic2(4,p3,states);
%             constraints=[constraints;temp2(:)'];
%             equals=[equals;1];
%         end
%     end
% end
% Tässä on apufunktio, joka tekee saman asian kuin ylläolevat sisäkkäiset
% loopit, mutta ottaa tilojen lukumäärän syötteenä, ja skaalautuu siten
% joustavammin
[Aeq, beq] = conditional_distributions(states);
constraints=[constraints;Aeq];
equals=[equals;beq];


%tässä sitten cross impact welhoutta
% prob(a)*Cab=prob(a|b)
% N:ssä dimensiossa täytyy ottaa huomion vanhempien yhteisjakauma, jotta
% saadaan laskettua ehdollinen todennäköisyys vain yhden vanhemman suhteen
%Cijkl*m(j)*p(kl)=cdist*(indexmagic(ijkl).*scenp(ijkl))
C=[];
d=[];
sceExt = ones(states(1),1)*scenp(:)';
for pi=1:nparents
    for sp=1:states(pi+1)
        psp = indexmagic2(pi,sp,states(2:end))'*skenprob2(:);
        for s=1:N
            C=[C;(indexmagic2(1,s,states).*indexmagic2(1+pi,sp,states).*sceExt(:))'];
            d=[d; C0(s,sum(states(1:pi-1))+sp)*m(s)*psp]; %#ok<*AGROW> 
        end
    end
end
C(end+1,:)=1E-3;
d(end+1)=0;
%C=sparse(C);
%constraints=sparse(constraints);
%ratkaistaan paras pns-sovitus
[x,resnorm,residual,exitflag,output,lambda]=lsqlin(C,d,[],[],constraints,equals,zeros(1,prod(states)),ones(1,prod(states)));
%sum(C*x)
%sum(d)
%residual
%x=reshape(x,states);
toc
end

%%% Apufunktioita
function [Aeq, beq] = conditional_distributions(states)
npstates=prod(states(2:end));
beq=ones(npstates,1);

Aeq=zeros(npstates,prod(states));
for i=1:npstates
    %iind=sum(states(1:i-1));
    Aeq(i,((i-1)*states(1)+1):(i*states(1)))=1;
end

end
