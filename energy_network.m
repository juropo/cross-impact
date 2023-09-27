% Tässä kuvaillaan Bayes-verkon perusrakenne, eli mistä solmuista on kaari
% mihin solmuihin, ja mitkä ovat solmujen nimet.
% CIindex kertoo, mistä kohtaa cross-impact matriisia mikin
% epävarmuustekijä löytyy, jos ne jostain syystä ovat eri järjestyksessä
nodes=struct("name","1. Nuclear power","parents",[],"children",[2 3],"CIindex",1);
nodes(2)=struct("name","2. Renewables","parents",[1],"children",[3],"CIindex",2);
nodes(3)=struct("name","3. Energy Storage","parents",[1 2],"children",[],"CIindex",3);

crossimpacts = readmatrix('energy_example.csv',"FileType","text");
crossimpacts = 2*log2(crossimpacts+crossimpacts');

%tilojen todennäköisyydet yhdessä vaakavektroissa samassa järjestyksessä kuin nodes
%tämän ei ole mitään järkeä olla vain yksi vektori
probs={[0.5 0.3 0.2], [0.3 0.2 0.3 0.2], [0.5 0.4 0.1]};
% states on vektori, joka kertoo, kuinka monta mahdollista tilaa jokaisella
% epävarmuustekijällä on. 
states=[3 4 3];

jointpd=[];

% Käydään läpi epävarmuustekijät yksi kerrallaan
for i=1:length(nodes)
    % Haetaan ristivaikutukset vanhempien kanssa
    CI = parentCIs(i,nodes, crossimpacts,states);
    % Lasketaan yhteisjakaumasta vain vanhempien yhteisjakauma
    pr = parentPD(i,nodes,jointpd,states);
    % Lasketaan i:nnen epävarmuustekijän todennäköisyysjakaumat
    % ehdollistetttuna vanhempien tiloille.
    nodes(i).cdist=ls_bayes(CI,probs{i},pr,states([i, nodes(i).parents]));
    % Lasketaan uusi yhteistodennäköisyysjakauma, jossa i on mukana, äsken
    % lasketun ehdollisen todennäköisyysjakauman perusteella
    jointpd = updatePD(i,nodes, jointpd, states, nodes(i).cdist);
end

% Alla olevalla voi tuottaa Bayes-verkon GeNIe-ohjelmistoon https://www.bayesfusion.com/genie/
% genie_parser(nodes,states,"tiedoston_nimi.xdsl");

function [CI] = parentCIs(ni,nodes, crossimpacts,states)
% Palauttaa ristivaikutusmatriisista sen osan, joka sisältään ni:n
% ja sen vanhempien väliset vaikutukset
CI=[];
i=nodes(ni).CIindex;
for p = nodes(ni).parents
    pi=nodes(p).CIindex;
    CI=[CI,crossimpacts(sum(states(1:i-1))+1:sum(states(1:i)),...
        sum(states(1:pi-1))+1:sum(states(1:pi)))]; %#ok<AGROW>
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
        if length(states)>1
            pr=reshape(sum(pr,i),states);
        else
            pr=sum(pr,i)';
        end
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