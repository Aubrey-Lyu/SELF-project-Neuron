# SELF-project-Neuron: Causal Evidence for the Processing of Bodily Self in the Anterior Precuneus
---

The repository contains the group-based ROI of self-hot electrodes and self-cold electrodes in the posteromedial cortex, and some key customized codes for generating the results of paper.


## Description of the data and file structure
- ROIs

The two nifti files are group-based ROIs from the study. With brain stimulation, we found a distinct area in the anterior precuneus that can induce robust subjective changes that are related to self-location displacement (sometimes self-dissociation). These stimulation sites are labelled "hot", while the surrounding posteromedial (PMC, or otherwise addressed as posterior cingulate cortex/precuneus) sites that did not induce such an effect are labelled "cold". The method of extracting the native coordinates of these stimulation sites is presented in the extractSbjNativeCoord.m. By transforming the individual brain to the standardized MNI space, we generated the presented ROIs by summarizing the spatial location of these stimulation sites (using 4mm-radius sphere). The ROIs are plotted in the provided figure ROI_vis.png, where the yellow color indicates hot sites. We hope these ROIs will be helpful for people who are interested in studying self-related processing, PMC heterogeneity and many other related topics. 

- fMRI/codes:

Among the *.m files, the fMRINetMapNative.m is a customized class function for doing all the fMRI-related processings, including preprocessing and seed-based FC at the native space as reported in the paper. The usage of this class can be found in the FCconndiff_selfvscold_example.m, which is the code that is used to conduct the seed-based FC of hot/cold stimulation sites and the contrast between them at the individual level. The group-level comparison (including normalizing the first-level stats images) is presented in the Elec_ContrElec_glm2.m. The results of this analysis was reported in the Figure 2 of the paper. 

MI_HCP_stanfordCohort_Comparison.m is the code for formally comparing the seed-based FC between the HCP cohort and our cohort. This part of result is presented in the supplementary figure S5.

- fMRI/FC_results:

HCP100Cohort and StanfordCohort respectively contains the group-level seed-based FC results from human connectom project open-access data (N=100, unrelated healthy young adults) and from our own cohort (N=5) where the hot/cold electrodes were identified individually. The subfolders with a suffix of "selfHot"/"selfCold"/"selfContr" respectively means the FC of the hot sites for the bodily self, the cold sites for that and the contrast between them. Among the nifti files, the xxx_0001.nii and xxx_0002.nii are associated with positive and negative contrasts, i.e. positive and negative connectivities, while the prefix of "beta", "con" and "spmT" indicates the image type (i.e., spmT is t-score images, beta is beta-coefficient images). When the suffix of an image is thr, it is a thresholded image based on the cluster-level significance (i.e. significant clusters), and a "**_bin.nii" is a binary image of the thresholded image (i.e. significant cluster mask). The resulting significance table is presented in the *.csv files. 

- CCEP

The CCEP_result_barplots.Rmd has the code for generating the barplot presented in the Figure 4 of the paper, and the prelim_CCEP_explore5_2.m has the code for generating the brain 3D plots of the figure as well as the inflow/outflow CCEP videos.

- metaTable
  This is a metatable documenting the electrode localization information for all the CCEP pairs. Explanation of variables: Self_Hot(1: hot, -1: cold); MNIout_coord_x (coordinate in the MNI305 space, which is compatible with fsaverage surfaces. "out" refers to the stimulation channel. "coord_1,2,3" refers to "x,y,z" axes respectively. The coordinates refer to the middle point of the bipolar contacts.); LEPTO_x_coord_x (coordinates in the native freesurfer surface space); Desikan_Killiany/DK/Destrieux refer to different brain atlases (automated output from freesurfer). JP_label (Josef's manual labelling based on individual high-res T1 scans); Native_Coord_x (native coordinates in real-world space).



## Sharing/Access information

Links to other publicly accessible locations of the data:
  * https://github.com/DianLyu577/SELF-project-Neuron/tree/main
  * https://doi.org/10.5061/dryad.w6m905qv0
  * https://github.com/Aubrey-Lyu/SELF-project-Neuron
  * Lyu, Dian; Parvizi, Josef (2023), “Spatial Localization of Anterior Precuneus for Bodily Self Validated with Brain Stimulation”, Mendeley Data, V1, doi: 10.17632/7pzyy4g6bx.1



## Code/Software
MATLAB >= 2019a
SPM12 
