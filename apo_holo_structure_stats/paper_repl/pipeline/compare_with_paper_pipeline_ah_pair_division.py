# ale chainy nesedi, oni tam maji ruzne
# jak vyresit? Bud spustit primo IsHolo na chainech
# A taky isoformu (ale to nejde) - slo by to, kdybych prepsal main chain po tom filter-structures
# takze ve finale to budeme delat na urovni chainu
# proste to nebrat jako strukturu, ale brat chain, a strukturu ignorovat

"""

1. vyzkousim, jak funguje is_holo, tak jak mam nastaveny, a isoform - jak budou vypadat ty groups (apo-holo snad jen 2, ale u holo-holo muzu vyzkouset cele grupy)
dve veci
1. parovani chainu/struktur
- brat jako chain (jako v paperu)
- kombinace pak vyzkouset bud vsechny, nebo nejak vzit nejlepsi match
2. co je sekvence - loop pres vsechny observed residues -- to uplne neni dobre, kdyz bychom unobserved primo ignorovali. Neco by tam mohlo být jiného.


Those proteins with little or no secondary structure have been removed.

Moreover, to avoid unresolved fragments of a structure and
misnumbering of residues, the amino acid sequences
used in this study were determined for the __fragments__ of  # takze tim chtej rict, ze to backbone-unobserved implicitne rozdelily na fragmenty? To by bylo dobře.
structures, for which at least backbone coordinates were
available.  - co si pamatuju, tak ja jsem to musel upravovat, protoze tam backbone nebyl
 - mozna, že měli LCS pak kratší, TODO zrušit to/porovnat když do LCS budu dávat jen tuhle sekvenci


To unobserved presne resi SIFTS:
The separate alignments for these segments were then merged together to assemble the complete alignment between the sequence of the observed
 residues from the PDB entry and the complete sequence of the protein used in the experiment.
 The latter sequence is shown in the SEQRES record in the PDB entry and does not have gaps reflecting unobserved residues.

Oni maj mapping observed --align-observed-segments-and-merge--> SEQRES -aligment-> uniprot

ja bych to mohl vyuzit, ze bych sel chain1 -> uniprot -> chain2. A musel bych kontrolovat, jestli se to mapuje na celej uniprot stretch (vcetne
unobserved) nebo ne
"""
