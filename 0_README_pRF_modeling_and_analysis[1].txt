pRF Modeling and Analysis 

These Jupyter notebooks perform the pRF modeling and analysis of the different sources of saliency (Eye Tracking First-Fixation and All-Fixations, and saliency estimations by DeepGaze IIE). 


*Information about the eye tracking experiment*

The eye tracking experiment involved 1,000 stimuli (retrieved from the Microsoft Coco Images Dataset) and 250 (1/4) of them were repeated. In the within-subjects design with randomized trials at a single time point, each participant completed 1,250 trials during data collection. The experimental conditions were the same as in the well-known experiment within the NSD dataset by Allen et al. (2022). Allen et al. (2022) systematically recorded only the voxel data of participants performing the same experiment while lying in an fMRI scanner. In contrast, our study captured the critical, missing ground truth eye-tracking data. We recorded eye movements from 23 participants as they viewed the same images used by Allen et al. (2022). The focus of this analysis is on the participants' eye movements, particularly fixations, during image perception.

Participants (*N*= 23) performed a continuous recognition task, which required them to remember and identify as many stimuli as possible. This indirect task design enables natural image exploration, yielding in more authentic eye-tracking data that better represents real-world visual perception.

------------------------- Notebook: 00_Creating_and_visualizing_roi_masks_04 -----------------------------------

This Jupyter notebook loads and visualizers the subject-specific ROI masks for the NSD subjects (n=8) and distinct ROIs.

Steps to run the code:

I. According to the selected subject number :  
II. It should load a list of the shared1000 Id  
III. Map the NSD image Ids to the trial Ids stored in nsd_stim_info_df_all dataframe
Prepartion of the pRF modeling 
1. Exploring the subject-specific angle data file (loaded from the NSD dataset), calculating the mid_vox 
2. Load the ROI masks from the NSD dataset 
3. Load pRF (eccentricity and angle) data and mask them with the subject-specific visual cortex masks
4. Save subject-specific ROI masks to disk
5. Visualizations of the relationship between the angle_data (polar angle) and the eccentricity. 
6. Function defintion load_prf_data 

------------------------- Notebook: 01_PRF_single_subj_shared1000_L78-iAll_Fix -----------------------------------

This Jupyter notebook perfroms the pRF modeling according to enterd NSD subject (n=8) and source of saliency (Eye Tracking First-Fixation and All-Fixations, and saliency estimations by DeepGaze IIE) 

Steps to run the code:

pRF Modeling of Beta Maps and Correlation of Modeled Data (Y hat) to Ground Truth Betas (Y)

Structure of this notebook:
- all_rois = contains all given rois 

*Part I*

1. Function definition and execution: “*get_saliency()*” (to load the saliency maps)
- 1.1- 1.3 Function call “get_saliency()” for one NSD image (for different sources of saliency predcitions (ground truth eye tracking data First-Fixation, All-Fixations, and aredictions according to the saliency model DeepGazeIIE))

2. Load betas (Y), all ROIs and pRF data from the NSD dataset
- 2.1 Function definition "*load_betas_concat(roi)*"
- 2.1.1 Function call *load_betas_concat()* for all ROIs
- 2.2 Function definition and execution *load_prf_data()*
- 2.3 Function definition and execution "*load_all_rois_prf_subject()*"
- *To load all ROIs from the corresponding pRF data with a function call of "load_prf_data()" and apply the selection criteria for valid voxels (r2_values> 0, only consider Voxels that have their pRF in the stim_sz_px + sigma/one standard deviation) to it.*
    - 2.3.1 Transfer **combined_prf_valid_df_subj1** DataFrame to a dictionary **prf_valid_dict_dfs-all_rois** for better accessibility

3. Computation of K via sequential approach and Einstein summation. Refine the valid voxels criteria to only pick those, where the PDF sums to one.

- 3.1 Load the pRF data for subject 1
- 3.2 Convert the polar coordinates of eccentricity data to cartesian coordinates of the shared1000 NSD images 
- 3.3 Generate pRF maps for each ROI -- Sequential approach: Slower way sequentially for each voxel 
- 3.4 Apply the Non_Zero_Mask for further selection of the valid voxels in df_valid to load only the beta data (Y) with the valid voxels in the correlation 
- 3.5 Visualization of examples of PDFs for each ROI
- 3.6 Compare K_mat of slower way sequentially processing for each voxel and the K resulting from faster computing via Einstein summation
- 3.7 Summing the probability densities across voxels (generated via Einstein summation) - Visualization

4. Function definition and execution: “*get_betas_from_numpy()*”

- 4.1 Function call ”*get_betas_from_numpy()*”
- *Load the three corresponding betas of one participant and calculate the mean beta map*
- 4.2 Load shared1000_list

5. Visualization of the relationship of sigma and eccentricity

6. Function definition and execution: “*match_salmap_and_betas_for_nsd_id()*”
    - Calls the functions *get_betas_from_numpy* and *get_saliency* at once, to match the outputs to the corresponding image
- 6.1 Function call “*match_salmap_and_betas_for_nsd_id()* for One NSD Id”
- 6.2 Function call “*match_salmap_and_betas_for_nsd_id()* for Muliple NSD Id's”
 

*Part II* Comparing pRF saliency based on ground-truth eye-tracking data to observed (fMRI betas) 

- A. Calculate Y_hat_and_get_betas and voxelwise correlation with fMRI betas
    - A.1 Function definition "*gen_corr_w_beta()* for at least two or more nsd ids stored in a list" 
    - A.2 Function call *gen_corr_w_betas()*" for all nsd_ids and v1
    - A.3 Function call generate correlation matrix and y hat for subject 1 and V4
    - A.4 Function call *gen_corr_w_betas()*" for all nsd_ids and V2

- B. Function call and execution *gen_corr_full()*" for all nsd_ids and all ROIs
    - B.1 Generate correlation corr_full for all ROIs/ a given ROI-list, for the given subject id and the given saliency model and the given nds_ids 
    - B.1.1 Function call: To generate correlations for all ROIs
    - B.2 Visualization of correlation
    - B.3 Generate correlation floc_faces
    - B.4 Generate correlation all_rois  


------------------------- Notebook: 02_Integration_of_DeepGaze_as_saliency_model_L10 -----------------------------------
This Jupyter notebook loads the shared1000 NSD images from the h5 file, initializes DeepGaze IIE and applies the saliency model to the shared1000 images. Returns image-specific saliency predictions according to DeepGaze IIE.

Steps to run the code:

Imports
Path Indications 
Load the TSV file of Shared1000 - NSD Images 
Installation/ Initialization of DeepGaze IIE 
Intialize DeepGaze IIE for running 
Preparation for correlation_w_fMRI: Load the image stimuli of shared1000 and load and mask the betas of the NSD dataset (resolution 1.8mm) for the relevant ROIs
Function defintion "gen_saliency()" to generate the model based saliency maps 
Now print applied DeepGazeIIE for all shared1000 images of the NSD Dataset 
Load the images from the h5py file 
Visualization 


------------------------- Notebook: 03_Analysis_corr_per_subj_L17 -----------------------------------
This Jupyter notebook analyses of the correlations of the modeled and ground truth fMRI responses within a single subject of the NSD.  

Steps to run the code:
- Initially: Set the sub_id and saliency model  
- Path Directions - NSD dataset or data directory.
- Sort the ROIs in the plot from early vision areas to higher vision areas 
- Load the necessary data and indicate the paths 
- Print the *p* values for the corrleation values of a subject (modeld beta maps to ground truth beta data)
- Apply Fisher's Z transformation before calculating the mean
- Calculate mean r using Fisher's z transformation and perform a one-sample t-test
- Plot the *r*-values for the significant rho values of df_roi_stats *(without FDR correction)
- Plot the ROI-specific correlation over all ROIs (with significance stars for significant ROIs)
- Apply FDR correction to *p*-values
- Comparison of the original p-values with the FDR-corrected p-values 

------------------------- Notebook: 04_Analysis_All_Subj_PRF_Model_Comparison_L20 -----------------------------------
This Jupyter notebook analyses of the correlations of the modeled and ground truth fMRI responses over all subjects of the NSD.  

Steps to run the code:

- Imports and path directions
- Sorting all ROIs according to their neuroanatomical order during visual saliency processing 
- Now mediate (arithmetic median) the results over the whole group. 
- One-Sample t-test for df_mean_roi_stats (DataFrame that stores over all ROIs, the corresponding Sal_model and the meaned r-values)
- One-Sample t-test for df_mean_roi_stats (DataFrame that stores for ROIs, the corresponding Sal_model and the meaned r-values)
- Visualize the r-values over all ROIs 
- Plot the mean correlation values for each ROI with significance stars and error bars
- Plot same barplot with confidence intervals and significance stars for ground truth eye tracking over all subjects.
- Plot same barplot with confidence intervals and significance stars for DeepGaze over all subjects. 
- Perform plots for single subjects 
- Rank df_mean_roi_stats for the highest and lowest ROI value 
- Paired t-tests: Is the correlation of the two saliency models (DeepGaze IIE or GroudTruthEyeTrackingData) significantly different?
- Repeated Measures ANOVA to evaluate the performance of all three models 
- Post-hoc: Paired t-tests (undirected), for each ROI, of the r values for DeepGaze and GroundTruthEyeTracking differ significantly across subjects.
- Model comparison averaged over all ROIs undirected paired t-tests 
- Now execution of directed t-tests
- Plotting of model difference (Mean r across all ROIs)

