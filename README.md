# the STOP g.APP

BETA release available here: https://stopgapp.elstondsouza.com/ <br />
This is a repository of scripts, documents and data associated with the "the STOP g.APP" RShiny app.
Variants that intersect stop codons have been strongly linked to a broad range of both rare, and common human genetic conditions, including cardiomyopathies, hearing loss and diabetes. However, interpreting the impact of any stop-loss variant is not straightforward as (1) the effect of addition of amino acids to the C-terminus on protein function is difficult to predict, and (2) if there is no alternative stop codon upstream of the polyadenylation signal the transcript may be degraded by non-stop decay. Differentiating variants that are likely to be deleterious from those that are likely tolerated often requires specialist analysis and skills as, at present, there are no tools available to assist in the interpretation of variants of this kind. 
<br />
The scale of this issue is illustrated in ClinVar where only 370 of 2,094 single nucleotide stop-loss variants (17.6%) have a confident classification (204 pathogenic and 166 benign).
<br />
To address this issue, we have developed ‘the STOP g.APP’, an R/RShiny web app that allows users to visualise potential effects of stop disrupting variants in any Ensembl v110 transcript. The simple, user-friendly interface enables users to investigate the impact of stop-loss variants (up to 50nt), providing information on downstream alternative stops, extension statistics, and whether non-stop decay is predicted. Importantly, the tool accounts for a change in frame when predicting the effect of insertion and deletion variants. The resultant annotations are provided as both a quick reference table, and an easily interpretable visualisation.
<br />
The STOP g.APP is a freely available open-source resource that can be used by both research and clinical communities to aid in the interpretation of stop-intersecting variants in human genes, and supplement interpretation using the ACMG-AMP guidelines. 
<br />
<br />
## Data: 
This folder contains all the background data used by the APP. Unfortunately due to issues with size we are not able to store the full data used to curate the app data on GitHub at this time. We do however provide all scripts used to generate these data in the ‘Data Scripts’ folder.
<br />
### Data Scripts: 
This folder contains all scripts used to generate the base data used by the STOP g.APP <br />
## WWW:
Logos and graphics used in the app
