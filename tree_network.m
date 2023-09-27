% Tässä kuvaillaan Bayes-verkon perusrakenne, eli mistä solmuista on kaari
% mihin solmuihin, ja mitkä ovat solmujen nimet.
% CIindex kertoo, mistä kohtaa cross-impact matriisia mikin
% epävarmuustekijä löytyy, jos ne jostain syystä ovat eri järjestyksessä

nodes=struct("name",int2str(1),"parents",[],"children",[],"CIindex",1);
for ci=1:81
    nodes(ci)=struct("name",int2str(ci),"parents",[],"children",[ceil(ci/3)+81],"CIindex",ci);
end
for i=1:27
    ci=81+i;
    nodes(ci)=struct("name",int2str(ci),"parents",find([nodes.children]==ci),"children",[ceil(i/3)+81+27],"CIindex",ci);
end
for i=1:9
    ci=81+27+i;
    nodes(ci)=struct("name",int2str(ci),"parents",find([nodes.children]==ci),"children",[ceil(i/3)+81+27+9],"CIindex",ci);
end
for i=1:3
    ci=81+27+9+i;
    nodes(ci)=struct("name",int2str(ci),"parents",find([nodes.children]==ci),"children",[ceil(i/3)+81+27+9+3],"CIindex",ci);
end
ci=81+27+9+3+1;
nodes(ci)=struct("name",int2str(ci),"parents",find([nodes.children]==ci),"children",[],"CIindex",ci);



%crossimpacts = readmatrix('cross-impacts-3d.csv',"FileType","text");
%crossimpacts = crossimpacts+crossimpacts';

%tilojen todennäköisyydet yhdessä vaakavektroissa samassa järjestyksessä kuin nodes
%tämän ei ole mitään järkeä olla vain yksi vektori
probs=rand(121,2)/2.01;
probs(:,3)=1-probs(:,1)-probs(:,2);

% states on vektori, joka kertoo, kuinka monta mahdollista tilaa jokaisella
% epävarmuustekijällä on. 3D-casessa se oli aina kolme, joten alla on vain
% vektori kolmosia. Koodi saattaa muuallakin olettaa, että tiloja on aina
% kolme, eikä vain katsoa states-vektorista, koska olen laiska
states=ones(1,length(nodes))*3;

jointpd=[];

% Käydään läpi epävarmuustekijät yksi kerrallaan
for i=1:length(nodes)
    % Haetaan ristivaikutukset vanhempien kanssa
    %CI = parentCIs(i,nodes, crossimpacts);
    CI = randi(7,[3,9])-4;
    % Lasketaan yhteisjakaumasta vain vanhempien yhteisjakauma
    %pr = parentPD(i,nodes,jointpd,states);
    pr = independentPD(i,nodes, probs, states);
    % Lasketaan i:nnen epävarmuustekijän todennäköisyysjakaumat
    % ehdollistetttuna vanhempien tiloille.
    % Tapa, jolla probs syötetään on vähän huono, koska se olettaa, että
    % jokaisella epävarmuustekijällä on kolme tilaa, mutta tuohon kohtaan
    % tarvitsee vain syöttää i:nnen epävarmuustekijän reunajakauma
    nodes(i).cdist=ls_bayes(CI,probs(i,:),pr,states([i, nodes(i).parents]));
    % Lasketaan uusi yhteistodennäköisyysjakauma, jossa i on mukana, äsken
    % lasketun ehdollisen todennäköisyysjakauman perusteella
    %jointpd = updatePD(i,nodes, jointpd, states, nodes(i).cdist);
end

% Alla olevalla voi tuottaa Bayes-verkon GeNIe-ohjelmistoon https://www.bayesfusion.com/genie/
% genie_parser(nodes,states,"tiedoston_nimi.xdsl");



function [CI] = parentCIs(ni,nodes, crossimpacts)
% Palauttaa ristivaikutusmatriisista sen osan, joka sisältään ni:n
% ja sen vanhempien väliset vaikutukset
CI=[];
i=nodes(ni).CIindex;
for p = nodes(ni).parents
    pi=nodes(p).CIindex;
    CI=[CI,crossimpacts(3*i-2:3*i,3*pi-2:3*pi)]; %#ok<AGROW>
end
end

function pr = independentPD(ni,nodes, probs, states)
% Laskee riippumattomien vanhempien yhteisjakauman kertomalla reunajakaumat
if isempty(nodes(ni).parents)
    pr=[];
    return
end
pr=1;
for pi=states(nodes(ni).parents)
    pr=pr(:)*probs(pi,:);
end
pr=reshape(pr,states(nodes(ni).parents));
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
% Lisätään ehdolliseen jakaumaan myös riippumattomat epävarmuustekijät
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
% Lasketaan yteisjakauma ehdollisen jakauman perusteella, ja järjestellään
% se oikein
skprob=ones(states(ni),1)*pdist(:)';
newScenarioP=(skprob.*reshape(cdist,states(ni),[]))'; %transpoosi lopussa heittää ensimmäisen epävarmuustekijän viimeiseksi
pr=reshape(newScenarioP,states([parents,ni]));

end