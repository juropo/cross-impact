%clear
diffsum=[];
jdiffsum=[];
for seed=1:10
rng(seed);
states=[3 3 3 3 3 3 3 3];
dist=random('exponential',1,states);
%dist=rand(states);
dist=dist/sum(dist,"all");

jpd=zeros(5,prod(states));

for i=1:length(states)
    nodes(i).name="variable"+i;
    nodes(i).parents=1:i-1;
    nodes(i).children=i+1:length(states)
    nodes(i).CIindex=i;
end

for iter=1:5

if iter==1
    perm=1:8;
    reverse=1:8;
else
    perm=randperm(8);
    [~,reverse]=sort(perm);
end

jointpd=[];
dist=permute(dist,perm);

for i=1:length(states)
    iind=sum(states(1:i-1));
    for k=1:states(i)
        probs(iind+k)=dist(:)'*indexmagic2(i,k,states);
    end
end

cmult=calculate_CI(dist,probs,states);
cmult=cmult+cmult';
crossimpacts=2*log2(cmult)'; 


for i=1:length(nodes)
    CI = parentCIs(i,nodes, crossimpacts);
    pr = parentPD(i,nodes,jointpd,states);
    [nodes(i).cdist, C0]=ls_bayes(CI,probs(1,3*i-2:3*i),pr,states([i, nodes(i).parents]));
    jointpd = updatePD(i,nodes, jointpd, states, nodes(i).cdist);
end

dist=permute(dist,reverse);
jpd(iter,:) = reshape(permute(jointpd,reverse),1,[]);

end

%tämän isomman ratkomiseen ei riitä muisti
jointpd2 = ls_joint(crossimpacts,probs,nodes,states);
jointpd2 = reshape(permute(reshape(jointpd2,states),reverse),1,[]);


%% 
%diffsum=[];
for i=1:4
    for j=(i+1):5
        diffsum(end+1)=sum(abs(jpd(i,:)-jpd(j,:)),"all"); %#ok<SAGROW>
    end
end

%jdiffsum=[];
for i=1:5
    jdiffsum(end+1)=sum(abs(reshape(jointpd2,1,[])-jpd(i,:)),"all"); %#ok<SAGROW>
end

end

function [CI] = parentCIs(ni,nodes, crossimpacts)
CI=[];
pr=[];
i=nodes(ni).CIindex;
foo=1;
for p = nodes(ni).parents
    pi=nodes(p).CIindex;
    CI=[CI,crossimpacts(3*i-2:3*i,3*pi-2:3*pi)]; %#ok<AGROW>
    %pr{foo}=probs(1,3*pi-2:3*pi); %#ok<AGROW>
    %foo=foo+1;
end
end

function pr = parentPD(ni,nodes, pdist, states)
% Laskee vanhempien yhteisjakauman summaamalla systeemin yhteisjakaumaa
% kaikkien ei-vanhempien yli
pr=pdist;
states=states(1:ni-1);
for i=ni-1:-1:1
    if ~any(i==nodes(ni).parents)
        states(i)=[];
        pr=reshape(sum(pr,i),states);
    end
end
end

function pr = updatePD(ni,nodes, pdist, states, cdist)
% Laskee päivitetyn yhteisjakauman uusimman lasketun ehdollisen jakauman
% perusteella
if ni==1
    pr = cdist;
    return
end
parents=nodes(ni).parents;
cdist=reshape(cdist,states([ni,parents]));
for i=1:ni-1
    if ~any(i==nodes(ni).parents)
        %mihin kohtaan puuttuva epävarmuustekijä kuuluu
        %huomaa, että ehdollisessa jakaumassa tarkasteltava muuttuja on
        %ensimmäisenä
        np=length(parents);
        perm=[1:i,np+2,i+1:np+1];
        %laske päivitetty ehdollinen jakauma riippumaton muuttuja ml.
        temp=cdist(:)*ones(1,states(i));
        %järjestele päivitetty ehdollinen jakauma
        cdist=permute(reshape(temp,states([ni,parents,i])),perm);
        parents=sort([parents,i]);
    end
end
skprob=ones(states(ni),1)*pdist(:)';
newScenarioP=(skprob.*reshape(cdist,states(ni),[]))'; %transpoosi lopussa heittää ensimmäisen epävarmuustekijän viimeiseksi
pr=reshape(newScenarioP,states([parents,ni]));

end

function ci = calculate_CI(dist,probs,states)
ci=zeros(sum(states));
for i=1:length(states)-1
    iind=sum(states(1:i-1));
    for j=i+1:length(states)
        jind=sum(states(1:j-1));
        for k=1:states(i)
            for l=1:states(j)
                ci(iind+k,jind+l)=(indexmagic2(i,k,states).*indexmagic2(j,l,states))'*dist(:)/...
                    (probs(iind+k)*probs(jind+l));

            end
        end
    end
end
end