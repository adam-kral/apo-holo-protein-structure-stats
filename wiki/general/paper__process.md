# Paper

## Filtering structures

[initial filter]
- X-ray, res <= 2.5 Å
- >50 residues (globally or one chain??)
- protein–DNA and protein–RNA complexes were excluded (interested in small molecules and short peptides)

[ligand filter], output: two groups: ligand-bound and ligand-free

- at least one ligand molecule is a ligand-bound, (all other structures are ligand-free)
    - >=6 non-H atoms (to exclude ions, water and very small molecules)
    - <= 15 standard or modified amino acids OK
    - single nucleotides OK

    - >= 6 residues in contact with the ligand as a whole (specificity)
        - observed as a result: each ligand on-hydrogen atom has six protein heavy atom neighbors (within a 4.5-Å cutoff) on average.

- those proteins with little or no secondary structure have been removed  (vždyť chtěli i sledovat změny v single domain proteins? Ale asi
ty se SS)

## Pairing structures

- make pairs of apo-holo at 100% sequence identity (that should allow one protein to be shorter than the other, however no alignment done
with gaps inside the sequence)
    - I would say there can be more pairs with same apo structure (pdb code), paired with different ligand-bound forms
        - as they say 'all possible pairs'

clustering procedure, using a cutoff of 35% sequence identity between clusters to remove redundancy
- The most representative pair from each cluster (closest to the cluster centroid) was selected to create the nonredundant set of ligand-
bound/ligand-free pairs used in this study (if more than one holo-structure for that sequence, choose a random one)

?nerelevantni Additionally those proteins (which? all of them?) were separately clustered by sequence similarity, and structure comparisons were done for the cluster centroids to investigate changes in a protein structure across different ligand-bound forms (tu ana
    - v results: dataset to nezminuji, kolik proteinu bylo vysledkem tohohle, pouze to 35 % identity a pak uz to delej na domeny
    - ani nikde dal v results jsem nevidel vysledky z tohohle

single domain proteins +
Multidomain proteins split into domains by Protein domain parser, result set:
- individual protein domains
- two domain arrangements - all possible domain pairs within a protein, (d1, d2), if their interface buried surface in both structural forms (ligand-bound and ligand-free) was larger than 200 Å

## Analysis

global structure similarity - RMSD of c_alpha atoms
local structure similarity - fraction of residues with same SS, as reported by DSSP (7+1 types, the +1 considered as random coil)

### motions of protein domains
- domain, principal axes (moment of inertia)
    - what exactly are the individual masses [36?]
        - [36] Thecenterofmassandthethreeprincipalaxes for domains definedforeachphosphorylasestructure were determined using only a-carbon coordinates fromtheindependentlyrefinedstructures (TableI).Eacha-carbonpositionwasequallyweightedinthemassdeterminations
            - c_alpha, equally weighted
    - isn't safer to do kabsch algorithm to make sure the axis correspond in corresponding domains?
        - for example - superimpose corresponding domains, average them, compute the principal axes and rotate them with the matrix
        - or maybe it was done? [36]?

- the single resulting angle is the largest difference (apo - holo) of angles between corresponding principal axes
    - largest |axis_i<dh1, dh2> - axis_i<da1, da2>|
    - jiny pristup -- maximalizace uhlu mezi u a Wu, pro u=u1,2,3 (nějaká báze U). Stačí minimalizovat SS <u, Wu> a získat tím u1, u2 pak s přidanou podmínkou <u1, u2> = 0 a u3 je pak jednoznačná, kolmá na obě a směr -- aby "pravotočivá soustava souřadnic"
        - dává pak podle mě větší smysl uvádět to jedno číslo, jako zástupné -- ne nějaký random rozdíl úhlů, ale úhel, kdy se okolo osy struktura otočila nejvíce (hinge osa?)
            - jinak se to myslim může rozložit, pokud určujeme úhly jenom mezi odpovídajícími vektory dvou random zvolených (hlavní osy vzhledem k momentu setrvačnosti) bázemi, a děláme ten rozdíl, nebo máme jen jednu zvolenou bázi, stejnou v d1 i d2, a zjištujeme úhly mezi bazickými vektory, když se báze ligand-pootočí (operátor W)
                - tam podle mě se stane to, že ten maximální z těch tří úhlů nebude stejný při volbě různých bazí (proto bych vzal tu, kde je), nebude správně reprezentovat otočení kolem dominantní hinge-osy a onen úhel, ale může se ten dominantní úhel "rozdrobit" do toho druhého, popř. i třetí hodnoty; nebude tak velký
        - tímto přístupem bych měl najít vlastně osu rotace (v podstatě jen jedna), to je ekv. nalezení vlastího vektoru k vl. č. 1 matice rotace. Ten úhel zjistim pak doplňím-li osu na ortonormální bázi a to zobrazení použiju na jeden z doplněných vektorů (jednu kterou) a pak třeba cos změřím úhel o kt. vektory zrotovaly

- SS similarity and RMSD of apo and holo, as well as comparing same protein with different ligands (different holo forms)


#### hinge or shear?
- The underlying assumption is that the hinge-bending motion of protein domains results in a considerable degree of bending as well as high RMSD between two conformers, whereas shear movements can be characterized by a detectable RMSD change and simultaneously, by relatively small angular movement of the domains’ principal axes.
jejich klasifikace (pozor, úhel je jejich - ten nejvetší rozdíl mezi úhly principal axis)
-  shear (RMSD > 1 and <2 Å)  =irrespective of bending (because in fact that would cause high RMSD and would by in hinge category)
-  hinge (RMSD > 2 Å, bending > 108°)

já bych ještě pro ověření shear typu na nějaké deformace samotné domény (to taky jde, viz paper FIgure 3 -- je zhruba stejné procento single domain protein i multi domain proteins s RMSD mezi 1-2, není důvod věřit, že všechny jsou shear type) zahrnul posun těžiště proteinu, jestli k tomu došlo. Popř. zkontrolovat podobnost posunových vektorů atomů. Nebo aspoň jestli RMSD ~= posun težiště. Další možnost -- čekal bych pak, že RMSD pak bude mít změněné

Vůbec mám trochu problém s (ne)analýzou toho Figure 2 (hlavně 0-2 RMSD)
- u multi-domain analýzy: pairů je 311, proteinů al 193. Tudíž se v grafu víc toto:  pokud se (jeden) ligand se naváže jen na jednu doménu (a není to large-scale movement), 2+ dalších zůstane netknutých a v grafu to bude vypadat (jak to vysvětlujou), jako že v multidoménových jsou domény víc netknuté a jde jen o large-scale movements (a s touhle "chybou" dál pracujou, že takový pohyby klasifikujou jako shear. Přitom to může být jen akomodace jedné domény, jako kdyby to byl single-domain protein).
- U RMSD 2+ to už vypadá spíš na to, co říkaj (domény multidomain proteinů jsou víc rigidní), ale, v závislosti na rozdělení počtu domén multidoménových proteinů, se individual domains mohou jevit v histogramu RMSD 2-3 výrazně nižší procento než single-domain-protein, ale musíme si uvědomit, že jsou to procenta PAIRS (a třeba u 4doménového proteinu by byla akomodující se doména jen ve 1/4 apo-holo pairů, tedy hypoteticky, kdyby nás zajímal počet proteinů, které akomodují (jeden) ligand jen jednou doménou, stejné zastoupení single-domain a 4domain proteins by se z pohledu PAIRS jevilo ve 4doménových čvrtinové!!!)  A ve vetšině případů ty proteiny mají v jedné struktuře ligand jeden (leda multimeric, ale ty, koukám v týhle studii nejsou, popř jen jejich monomery).

- ledaže by to byly ligand-bound ligan-free pairs z pohledu domén, ne proteinů (individual protein domain)


Large-scale domain movements kapitola
- filtrovali to teda (aby věděli, že je large-scale a ne jenom v jedný doméně třeba), nebo je tam ten problém jak popisuju výš


klasifikace sedí s -- Classification of macromolecular motions proposed by Gerstein and coworkers.15,21
- neni to to stejny treba? Podivat [15], [21]

http://molmovdb.mbb.yale.edu/molmovdb/

#### surface area vs motion types:

- The largest average bending and average RMSD was found for the relatively small interface surface (between 1000 and 1500 Å^2)
- 1500–2500 Å^2 both main types of motions
- More extensive interactions between domains
(2500–3500 Å^2 ) account for the total preclusion of
hinge-bending movements; however shear movements
still occur.
- no large-scale movements were
observed if the interdomain surface area was greater than
3550 Å^2 .

Those results suggest that some preferences toward the anticipate type of motion are deducible by the examination of a ligand-free structural form alone.


## Moje otázky

Jak funguje PDP (protein domain parser)
Molmovdb a klasifikace pohybů, jak např. definuje shear? Líp, než v tomhle paperu?
- dále pak popisujou asi shear: depending on whether or not they involve a parallel-plane sliding of one domain relative to the other. (ale sami to nepočítaj)

Které struktury uvažovali?

Jak spolu napárovali struktury?
- 100% sekvenční identity?
- musí mít i stejnou délku?
- jak se budou lišit výsledky – počet párů – když povolím buď různé délky nebo pár mismatchů?


Jsou v pdbe-kb uniprot-skupinách všechny pdb struktury?
    22.11.2020 uniprot_segments_observed 160059
    23.11.2020 pdb #structures 171313

    párovány nástrojem SIFTS  -structure-based annotations for proteins
    - original article https://academic.oup.com/nar/article/41/D1/D483/1072377
    - Structure Integration with Function, Taxonomy and Sequences
        -  semi-automated process to identify sequence cross-references from UniProtKB to the protein sequences in the PDB, and a fully automated process to generate residue-level mappings between the two sequences

        (i) mapping to a UniProt canonical protein sequence, unchanged compared to the previous implementation,
        (ii) mapping to all alternative isoforms of the canonical sequence
        (iii) mapping to sequences in UniRef90 clusters.

    original SIFTS
    - They must have a high level of sequence identity (ideally 100% but not below 90%);
    - The source organism must be identical or must have a common ancestor within one or two levels up to species level in the taxonomy tree.

    - **The process results in automatic identification of the correct UniProtKB cross-reference for 80–90% of the PDB entries.**
    - Sequences of biological origin not contained in the UniProtKB database are flagged for inclusion into UniProtKB.


    V uniprotu v kategorii existence:"evidence at protein level" 276K,
        of that 171K unreviewed – no publication, just from large proteomic data for example),
        105K reviewed – have a publication declaring that the protein exists (that would be all proteins with a determined structures I would say)

        -> <105K uniprot groups

Jak vypadá uniprot skupina, mají všechny její struktury 100% sekv. identitu i délku?
- Ne, ale možná, že v praxi často ano.


[TODO] ještě prohlédnout další papery, (z těch jiných db? jestli tam dělaj podobnou analýzu)