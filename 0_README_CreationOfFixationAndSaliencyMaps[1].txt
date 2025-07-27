Creation of Fixation and Saliency Maps for Visual Saliency Analysis of ET COCO Experiment

*Information about the eye tracking experiment*

The eye tracking experiment involved 1,000 stimuli (retrieved from the Microsoft Coco Images Dataset) and 250 (1/4) of them were repeated. In the within-subjects design with randomized trials at a single time point, each participant completed 1,250 trials during data collection. The experimental conditions were the same as in the well-known experiment within the NSD dataset by Allen et al. (2022). Allen et al. (2022) systematically recorded only the voxel data of participants performing the same experiment while lying in an fMRI scanner. In contrast, our study captured the critical, missing ground truth eye-tracking data. We recorded eye movements from 23 participants as they viewed the same images used by Allen et al. (2022). The focus of this analysis is on the participants' eye movements, particularly fixations, during image perception.

Participants (*N*= 23) performed a continuous recognition task, which required them to remember and identify as many stimuli as possible. This indirect task design enables natural image exploration, yielding in more authentic eye-tracking data that better represents real-world visual perception.

------------------------- Notebook: CreationOfFixationAndSaliencyMaps_26 -----------------------------------

This Jupyter notebook Creates fixation and saliency maps based on the eye tracking data collected in the mentioned experiment.

Steps to run the code:

Run Notebook cells linearly one after the other, because they are dependent on the previous cells.

I. Import the necessary libraries
II. Base directory, data directory and save directory setup

1. Complete run through all 1000 images to create fixation and saliency maps  
2. Further cells of this notebook are useful, if you want to (re)create the visualization of a certain image: 

- Fixation Map of all Fixations (Option: Enumerate fixations)
- Fixation Map of only the first fixations during the first presentation of the image for the participants (Option: Enumerate fixations)
- Saliency Map of all Fixations 
- Saliency Map of only the first fixations during the first presentation of the image for the participants
- Print the fixations on a scatter plot, representing the dimensions of the experimental screen 



------------------------- Notebook: FilterDataFrame_Eliminate_Fixations_50ms -----------------------------------


This Jupyter notebook includes an option to eliminate the rows with fixation_durations < 50ms (n= 7,312 Fixations) and therefore filter the DataFrame as suggested by Rothkegel et al. (2019). 

Steps to run the code: 

Run Notebook cells linearly one after the other, because they are dependent on the previous cells.

I. Import the necessary libraries
II. Base directory, data directory and save directory setup
III. Test accessibility of concatenated DataFrame 'df_expanded_data_subjs_concatenated_all.csv' in the save_dir

 
A. Eliminate the rows with the fixation_durations <50ms (.05ms) from the DataFrame
