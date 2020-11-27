# Input: creating groups of structures, usable APIs

## sifts mapping
### Uniprot sequences intro 

Uniprot entry
- "Each entry corresponds to a single contiguous sequence as contributed to the bank or reported in the literature." -- https://www.uniprot.org/docs/userman.htm

Sekvence v uniprotu mohou mít různé identifikátory:

primary accession, uniprotkb entry  ( P04150)
isoform accession  (P04150-2, primary+dash+number)
chain (sequence after PTM) identifier

určují hierarchické entity, od otce: (primary, isoform, chain)


ALE isoform teoreticky může mít víc otců:
#### Isoformy 
When alternative sequences _differ significantly_, we create separate entries for them and list all isoforms in each of them.
- isoforms produced from a single gene listed in one entry may have identifiers derived from different primary accession numbers.
- differ significantly -- takže mě to nezajímá, proteiny od dvou takových sekvence bych stejně nesuperimposoval 

Co je isoforma? -- https://www.uniprot.org/help/alternative_products
- alternative protein sequences (isoforms) that can be generated __from the same gene__

Jak vzniká?
- a single or by the combination of up to four biological events (alternative promoter usage, alternative splicing, alternative initiation and ribosomal frameshifting)



SIFTS api dokážou vrátit i chain identifier.

Mohou se izoformy lišit sekvencí na stejných místech v genomu? Snad ne. Pak by se vytvořil nový entry.

Důležité je pro mě, zda stačí využít ty primary accession, jelikož pouze ty dostanu v těch batch csv..
- jestli ty observed segments jsou 



### K čemu SIFTS mapping?
To jsou ty struktury, které mě budou zajímat + (snad) tím získám mapping  mezi apo a holo strukturama stejnejch proteinů (chainů)
- no mapping tim získám možná jen na úrovni 

### Jak vzniká?
- https://www.ebi.ac.uk/pdbe/docs/sifts/methodology.html
- depending on _taxonomy_ + _sequence identity_ 
    - (they cleaned up+structured the taxonomy data for the pdb entries)
- Sequences of a structure in the PDBe may represent either the native protein sequence or that of an engineered mutant or other variant, during the automatic procedure the criterion for assessing sequence identity was that there should be 95% or higher agreement between the sequence of a protein structure and the corresponding sequence in UniProt.
- mapping taxonomy rule
    - allows the taxonomy ID for the two entries, PDBe and UniProt, to be the same or to have a common parent within one or two levels up the taxonomic tree (since structure is more conserved than sequence during evolution time)
        - but I guess it's ok, as the sequence still has to be 95 % identical and the best uniprot match is chosen anyway
        
- (so, as both conditions have to be met, they probably chose the second one, in the process of computation of te mapping, as it's computationally less expensive)
### 
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
        
        
 ### files
ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/
ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz

Mělo by obsahovat 80-90 % pdb entries (2012)
- 95% sekv. ident. (90% když není nalezen match a má to být doděláno manuálně)
- často totiž mají proteiny engineered mutace (to asi neni v uniprotKB jako sekvence...)

#### files content
Co znamená None v PDB_BEG PDB_END? to me asi nezajima, zajima me jen RES_BEG a RES_END... (to je SEQRES a to PDB je pdb sequence, ale to je ta author sequence, ktera nikoho nezajima?)
jak se liší pdb_chain_uniprot (pdb to uniprot) a uniprot_segments_observed (uniprot observed in pdb) (větší)
- obojí má někde rozdělený úsek na fragmenty, kde v tom druhém není rozdělený
- PCU má None v PDB_BEG/PDB_END, USO nemá
- ani jedno nemá uniprot_id isoformy, takže to nestačí (isoformy se vždy liší?) a musím dělat dotaz na api, které je má (jednotlivě)
- to ftp:xml (jednotlivé struktury) csv etc. je asi jenom to stare api, ktere nevypisuje ty isoformy?
- možná když použiju to api s isoformami na mapping residues, nebude mě zajímat, jak se tyhle dva fily liší


- pozor -- možná uniprot_segments_observed nepočítá s isoformami a pdb_chain_uniprot jo, i když tam nemá přímo id ty isoformy, ale jen primary accession
    - nevim, ale v pdb_chain_uniprot jsou i unobserved (nemají ATOM recordy)
pdb_chain_uniprot
- čím je pokryto SEQRES (tzn. je-li v seqres něco navíc oproti uniprotu, bude tam víc záznamů rozdělených)
- rozdělené -- neboli ve struktuře je něco navíc
uniprot_segments_observed
- rozdělené -- v uniprot je něco navíc (unobserved residues as ATOM coordinates in structure)
- pokud je záznam stejný, je to kandidát pro moji skupinu. Pokud je jiný, automaticky nebude mít 100% identitu a není kandidátem
    - případně pro fuzzy by šlo groupovat jen podle UniprotKB primary accession a skupiny tvořit sám
    - ale není to jistý, musím porovnat všechny residues (co když nejsou v některé struktuře všechny zachycené, ale asi by tam měly být -- dle author_seq), co když tam neměly být?

3k23 -- divný,

__vedlejší pozn.__
jak vzniká entity_poly.pdbx_seq_one_letter_code? Ta sekvence v mmcif? Je tam exp. základ, jen není možné získat přesnější coordinates, nebo to není vůbec pozorováno?
vs auth_seq_id (asi myšleno do něčeho jako uniprotKB sekvence)


[todo] hledám co měli za uniprot sekvence (entrys consensual?)


### API services

https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/:accession
- získat nejlepší struktury s dobrým coverage a resolution (co když je mezi unp_start a unp_end start a end SNP??)
- ale je to po jednom
- to ftp:xml csv etc. je asi jenom to stare api, ktere nevypisuje ty isoformy?

https://search.rcsb.org/index.html#search-attributes
any mmcif dictionary value (and additions rcsb)
- has specifiable return type (entity_type)



PDBe apis
- nikdy celé struktury (nenašel jsem), jenom metadata a to ne všechny najednou, jaké chci.
- entry-based api
- aggregated api -- metadata o víc strukturách (např nejlepší struktury pro uniprot id)


## Creating a structural model

Biopython
- lze nacist bez headeru, o neco rychleji (cte pouze atom lines) a header jen par polí, málo pro generic filtering
- pokud by nás zajímal header, lze nacist mmcif do dictu. Pak by se ten dict dal vyuzit pro build struktury, ale nemaj to tak udelany

Biopython rychlost
- sám biopython.PDB uvádí 1.5 s/ structure (800 proteins each from a unique SCOP superfamily)
    - tzn. 80 hodin pro 200K struktur, a to ještě žádná analýza (aktuálně: vypadá to, že min. 10 h pro 100K struktur podle toho níže)
    -
- https://github.com/project-gemmi/mmcif-benchmark#creating-a-structural-model
- např. gemmi-structure 16x rychlejší na struktuře s normální velikostí 836 kB (75 percentil)  (biopython 0.5 s)
    - 30x rychlejší na obří struktuře (HIV-1 capsid 250 MB)  (biopython 128 s)
    - ale developed primarily for use in macromolecular crystallography
- nebo ihm (integrative-hybrid modeling) https://github.com/ihmwg/python-ihm, podporuje i BinaryCIF
    - co umí?
    - c-extension parsování
    - mmcif into a hierarchy of classes




Vstup do dalšího algoritmu bude:

struktury, mapping struktur apo-holo (nebo skupiny k uniref jako tady)

ty namapované ještě (už s jejich pdb fily) zkontrolovat, jestli mají 100% sekv identitu (jako v paperu).



Pokud možno a výhodno offloadovat část filtrování na server.

ty struktury dál vytřídim (jak? Na počítači, nebo dotazem do server databáze) na ty dostatečně kvalitní (resolution), dostatečně residues,  a s ligandy (dostatečně hetatm nebo peptid = krátkej AA chain??)


clustering si pak dělám sám (používám sekvence z pdb asi), nebo jsou nějaký auto-verifikovaný z pdb (auto-annotation)??


## SIFTS download and processing
### Prvotní analýza (počet, velikost skupin)
skupiny jsem dělal tak, že určitý pdb chain se s jiným musel naprosto shodovat v jejich uniprot mappings (observed, nebo i klidně merged z toho druhého souboru). Nevím ale, jak přesně ty observed fragmenty vznikly (asi pořád platí, že tam stačí celkově 95% identita, toho merge fragmentů, takže AK nejsou nutně stejné -- teoreticky by konce těch fragmentů mohly být fuzzy a možná i u struktur, které bych chtěl v páru porovnávat by se trochu lišily (ale chci porovnávat 100% sekv. identitu, tak by fragmenty asi měly být stejné -- vznikly za stejných podmínek-- stejná vstupní sekvence))
chci říct, že pro tu analýzu možná budu chtít groupovat jen podle uniprotkb_id a ne nutně podle těch fragmentů/merged mappings. Pak ale ještě budu kontrolovat opravdu sekv. identitu při načítání párů struktur /to stejně/.
[todo] .describe dataset pouze grouplý přes uniprotkb_id

 