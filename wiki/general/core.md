
# Architektura poznámky

jednotlivé analýzy mohou na sobě záviset ( dependency graph). Taky by se to mohlo ignorovat a pokaždý počítat zvlášť. Ale: co přístup, kde se budou ty analýzy cachovat.

A pozor, budou mít stejnou délku i 100% identitu stejné uniprot_id? Pak by bylo matchování residuí (domén atp.) triviální

Ale co analýzy, které je možné provést různými způsoby (prográmky), např. detekce domén. Tam by se měl cachovat výsledek. Ale co když v nějaké analýze bych chtěl domény z jiného prográmku?
get_domain. Zrovna domény ale nejsou na té zajímavé úrovni.

RMSD -- ty získám pro jednotlivé domény (úroveň analýz -- jednotlivé domény). Když ale budu analyzovat páry domén, budou mi stačit výsledky z těch jednotlivých: rotační matrix první domény (R1) a druhý (R2). Obojí získané z nezměněných struktur. Tak zarovnání páru na tu první -- otočení R1. Pak je neznámá R otočení druhé domény tak, že R2 = R1*R; R = R1^(-1)*R2

[Oni ty domény ještě snad averagovali (to mělo možná smysl u principal axes metody, ale tady by to nic neměnilo -- musel bych je otočit zas zpět, ze zarovnané polohy a úhel by byl tedy stejný, šlo by o stejné otočení))]

Taky určit, co bude output a co ne (může mě zajímat jen RMSD, ale budu cachovat všechno, i rotation matrix).

Úroveň analýz páry domén -- bude se mi hodit rotační matrix z těch jednotlivejch domén. (Možná to předám rovnou do workeru týhle úrovně.) Protože analýza těch párů je založená na správné vzájemné orientaci. (Vlastně ne -- např. contact area nepotřebují otáčet., ani třeba SASA v místě ligandu)

Filter může obsahovat analýzu!!!
-- interface buried surface --> min. 200 Å^2 pro analýzy large-scale domain movements
-- a je blbost počítat analýzu znova, byla-li již vypočítaná ve filtru



Vždy pár apo-holo, popř. holo-holo. Říkejme tomu pdb_struct (ne izoforma).
[todo] Co když má víc Modelů struktura? Může se stát, že jeden bude třeba apo a druhý holo?
RMSD v rámci jedné struktury nedává smysl (ledaže by to byl multimer -- otázka jestli analýza samotných multimerů je něco, co má tento program podporovat -- asi to půjde udělat, ale nebude to v pipelině. A není to pak relevantní (nacachované data) pro 2struct pipelinu)

- je nějaká featura relevantní z 1struct analýzy pro 2struct analýzy (kde se jedna struktura může opakovat několikrát např. v apo-holo, holo-holo párech)?
    - centroid domén (pak pro počítání RMSD) (v různých struct různý)
    - sekundární struktura pro jednotlivé residues (pak pro počítaní % identity SS)


Multimery
-  jak napárovat odpovídající chainy pokud homomultimer? -- Mohlo by to být jedno, pokud má multimer určitou symetrii. Obecně to jedno asi neni. Možná pomocí zarovnání? Když se zarovnaj na sebe.
- jinak je analýza stejná jako výše, až na to, že za domény dám naopak chainy

- pro každý multimer se pak budou analyzovat jeho chainy zvlášt na doménovou analýzu. Ale, odpovídající chainy u homomultimeru asi zas všechny možnosti zkusit

- u multimerů prostě musí být vrstva, která, v případě výskytu stejných monomerů rozhodne, které chainy dvou struktur sobě odpovídají/zkusí všechny možnosti (můj default třeba).
    - tohle (párování) pak půjde pro analýzu vzájemného postavení podjednotek i pro analýzu jednotlivých podjednotek (domén v nich).


Analýzy, (kt. jsou zamýšlené pro výstup) musí být serializovatelné. Pak se ještě musí připojit metadata úrovně zpracovaní -- o jakou strukturu a chain jde (popř. kombinace).

Ideálně bych chtěl udělat RDF výstup. Ale nenašel jsem nějakou Code-generation knihovnu podle owlu, která by zajistila, že jsem v kodu nic neopomněl a umožnila mi využít owl/shacl shape file jako popis výstupu programu. Je ovšem pravda, že výstup se může měnit podle nastavené pipeliny. Tak aspoň owl, aby se to dalo propojit, s tim, co už definice má (chain, structure, rmsd atp.), abych to nemusel popisovat já.



BinaryCIF
- neni pro pdbe data
- parsovani neni o tolik rychlejsi (ale js, ne python, ten ma c-extension)
- je asi 3-4x mensi, pri komprimaci deflate 5x nez text cif


# Pipeline

pipeline naše konfigurace

workflow

get structures with same uniprot chains
filter them (quality etc.)
classify apo/holo
make pairs (apo-holo, holo-holo)
    feed them to bistructure analysis
    identify domains (in one or both, if 100% identity, one is OK)
        bistructure-bifragment analysis



analyzers - filter uses them, a classifier is one, bistructure analysis, bistructure-bifragment analysis
    should be somehow cached (probably only when processing th
feeders
    take list of structures (potentially classified etc.)
    run analysis with structures (pairs)

cache is always for a pdb structure (and it's fragments). Will be a structure used in analysis multiple times?
    yes, in the uniprot-chain_id-class, only the chain
    if I wanted also multimer (as whole) analysis, I would have to detect which structures are of the same multimer (composed
of same chains might be enough), let's say we have the mapping. Than, again only within the same-bio-multimer group one structure could
be used multiple times.
    ok, finally if we combined those analyses, could they share some cached data (be sharing beneficial?)
        probably not, because those would not be the same structures (above only monomeric proteins)
        but we could do the monomeric analysis for (corresponding) monomers
            then, it could be cached (centroid of monomers, secondary structure of domains -> monomers)
                but seems not that important..., again optional, just see, if it would be simple to do (esp. if the caching would work
                for the monomeric analysis
            yes, the multimeric analysis would be prob. same as monomeric, but in multimers it would add the higher-scale data layer


        the caches could be valid as long as the structure object (loaded struct) exists -- (one would in the feeders maximize the life 
of the structure,  not load it again for no reason


analyzer
- the result should be serializable (or could be on result class)
- it should be cached
-

analyzers may depend on each other
- would there be a