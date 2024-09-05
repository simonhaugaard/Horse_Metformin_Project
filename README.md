# Bulk RNA-seq analysis of horse AF data

In this study, we induced atrial fibrillation (AF) in a cohort of twenty retired racehorses and followed them over 4 months. Half of them were treated with metformin, other half placebo treated. We also included a sham group comprising four healthy horses in normal sinus rhythm. Finally, there were six horses with naturally occuring AF presenting at the Large Animal Teaching Group for treatment, from which we also had biopsies from both atria. 
- **Placebo group**: Four months of AF, placebo treated (n = 10)
- **Metformin group**: Four months of AF, metformin treated (n = 10)
- **Sham group**: Sham-operated in normal sinus rhythm, no treatment (n = 4)

After 4 months, we acquired biopsies from the left and right atrium of the heart. Here, we wanted to assess differences in atrial gene expression between the groups.

Uniquely to proteomics, we also had baseline biopsies from right atrium before any interventions allowing for a "before/after" ("Timecourse") comparison. 

Our primary objective is to uncover how four months of AF alters the atrial transcriptome, and assess the potential effects of metformin in ameliorating such changes. 

# Directory Description
RNA-seq: (Short read bulk sequencing) This folder contains analysis related to the DGE-analysis comparing  all three groups. 
- **data:** Count matrix, metadata, gene annotation
- **analysis:** Differential gene expression analysis followed by gene enrichment

Proteomics: (LC/MS) This folder contains analysis related to the DGE-analysis comparing  all three groups. Directory as for RNA-seq, but contains to main branches 
- **Terminal:** Data & Analysis from terminal (4month) biopsies, comparable to RNA-seq
- **Timecourse:** Data & Analysis from terminal (4month) biopsies and baseline (0months) biopsies - here we look at the changes over time. 


# LC/MS Proteomics analysis of horse AF data
Here, we performed LC/MS proteomics on the first three groups (not natural AF). Other than comparing the metformin_vs_control & placebo_vs_sham, we also had samples transvenously obtained baseline samples from before AF-induction. 
This allowed for a timecourse comparison... 
