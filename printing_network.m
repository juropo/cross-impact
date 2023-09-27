% Tässä kuvaillaan Bayes-verkon perusrakenne, eli mistä solmuista on kaari
% mihin solmuihin, ja mitkä ovat solmujen nimet.
% CIindex kertoo, mistä kohtaa cross-impact matriisia mikin
% epävarmuustekijä löytyy, jos ne jostain syystä ovat eri järjestyksessä
nodes=struct("name","1. Global industry growth","parents",[],"children",[2 3 4 5 7 8],"CIindex",1);
nodes(2)=struct("name","2. Printing speed compared to present","parents",[1],"children",[3 4 11],"CIindex",4);
nodes(3)=struct("name","3. Printing cost compared to present","parents",[1 2],"children",[4 8 9 10],"CIindex",3);
nodes(4)=struct("name","4. Finnish industry growth","parents",[1 2 3],"children",[5 6 7 8 9 10 11],"CIindex",2);
nodes(5)=struct("name","5. Number of graduates with 3D-printing expertise","parents",[1 4],"children",[9],"CIindex",5);
nodes(6)=struct("name","6. Legal regulation of 3D-printing in Finland","parents",[1 4],"children",[8 9 10],"CIindex",6);
nodes(7)=struct("name","7. Standardization of processes and models","parents",[1 4],"children",[9 11],"CIindex",7);
nodes(8)=struct("name","8. Use of 3D-printed objects in FDF","parents",[3 4 6],"children",[9 10 11],"CIindex",8);
nodes(9)=struct("name","9. FDF access to 3D-printing model files","parents",[3 4 5 6 7 8],"children",[10 11],"CIindex",9);
nodes(10)=struct("name","10. FDF 3D-printing spare parts in peacetime","parents",[3 4 6 8 9],"children",[11],"CIindex",10);
nodes(11)=struct("name","11. FDF 3D-printing spare parts in crisis times","parents",[2 4 7 8 9 10],"children",[],"CIindex",11);

crossimpacts = readmatrix('cross-impacts-3d.csv',"FileType","text");
crossimpacts = crossimpacts+crossimpacts';

%tilojen todennäköisyydet yhdessä vaakavektroissa samassa järjestyksessä kuin nodes
%tämän ei ole mitään järkeä olla vain yksi vektori
probs=[0.3 0.5 0.2, 0.4 0.5 0.1, 0.5 0.4 0.1, 0.3 0.5 0.2, 0.2 0.6 0.2, 0.05 0.9 0.05, 0.35 0.45 0.2, 0.1 0.5 0.4, 0.2 0.7 0.1, 0.7 0.29 0.01, 0.45 0.45 0.1];
% states on vektori, joka kertoo, kuinka monta mahdollista tilaa jokaisella
% epävarmuustekijällä on. 3D-casessa se oli aina kolme, joten alla on vain
% vektori kolmosia. Koodi saattaa muuallakin olettaa, että tiloja on aina
% kolme, eikä vain katsoa states-vektorista, koska olen laiska
states=ones(1,length(nodes))*3;

jointpd=[];

% Käydään läpi epävarmuustekijät yksi kerrallaan
for i=1:length(nodes)
    % Haetaan ristivaikutukset vanhempien kanssa
    CI = parentCIs(i,nodes, crossimpacts);
    % Lasketaan yhteisjakaumasta vain vanhempien yhteisjakauma
    pr = parentPD(i,nodes,jointpd,states);
    % Lasketaan i:nnen epävarmuustekijän todennäköisyysjakaumat
    % ehdollistetttuna vanhempien tiloille.
    % Tapa, jolla probs syötetään on vähän huono, koska se olettaa, että
    % jokaisella epävarmuustekijällä on kolme tilaa, mutta tuohon kohtaan
    % tarvitsee vain syöttää i:nnen epävarmuustekijän reunajakauma
    nodes(i).cdist=ls_bayes(CI,probs(1,3*i-2:3*i),pr,states([i, nodes(i).parents]));
    % Lasketaan uusi yhteistodennäköisyysjakauma, jossa i on mukana, äsken
    % lasketun ehdollisen todennäköisyysjakauman perusteella
    jointpd = updatePD(i,nodes, jointpd, states, nodes(i).cdist);
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