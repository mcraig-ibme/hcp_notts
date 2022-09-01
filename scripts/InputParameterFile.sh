# The purpose of this file is to give location of the subjects' folders and the names of the images. Also, the scanner parameters and distortion correction methods are set up in this file as well as the location of the templates (which normally do not need to be changed). 

# Change the scan parameters below to match your acquisition, these are set to match the HCP Protocol by default

#If using gradient distortion correction, use the coefficents from your scanner
#The HCP gradient distortion coefficents are only available through Siemens
#Gradient distortion in standard scanners (even 7T scanners) is much less than for the specific HCP Skyra, and this correction step can typically be skipped.


### FOLDER LOCATIONS ###

# To use this pipeline the folders and images must be arranged as follows:
#   - a single top-level directory must contain one subdirectory for each subject
#   - within each subject directory there must be all the necessary image files
#   - image files cannot be within separate subdirectories
# "StudyFolder" describes the single top-level directory containing the individual subject directories
# "Subjlist" contains a list of the subject directory names
# Example of folder structure: 
#     /home/fs0/myname/scratch/mystud/CON077/
#     /home/fs0/myname/scratch/mystud/CON078/
#  would be described by:
#   StudyFolder=/home/fs0/myname/scratch/mystud/
#   Subjlist="CON077 CON078"

StudyFolder="/gpfs01/share/HCP/HCPyoung_adult/Structural" 

# Write space delimited list of subject IDs. If the filename of images contains the subject ID (from the Subjlist) use the wildcard %subjectID% instead of writing the actual subject ID, e.g. if the filename on disc is "T1_CON077.nii.gz" use the form "T1_%subjectID%.nii.gz" in the lines for Input images in this file
Subjlist="100206"  	
  

##############################
###### Structural input ######

#Naming conventions:	StudyFolder/SubjectFromSubjlist/image.nii.gz 
 
  #Define inputs as they are in the study, replace subject ID number (if is part of the filename) with a wildcard %subjectID%
  # T1w
  #T1wInputImages="T1A.nii.gz@T1B.nii.gz"
  T1wInputImages="T1w.nii.gz"
  
  # T2w (if it exists, otherwise set T2wInputImage to "NONE") 
  #T2wInputImages="T2A.nii.gz@T2B.nii.gz"  # filename of T2w image or "NONE" if not used 
  T2wInputImages="NONE"

  # B0 Distortion correction of structural
  AvgrdcMethod="NONE" #Averaging and readout distortion correction methods: "NONE" = average any repeats with no readout correction "FIELDMAP" = average any repeats and use field map for readout correction "TOPUP" = average and distortion correct at the same time with topup/applytopup only works for 2 images currently

       ### If you use AvgrdcMethod="FIELDMAP" change variables in this part. Otherwise ignore. Fieldmap can be used only for Siemens scanners.
       FmapMagnitudeInputName="FieldMap_Magnitude.nii.gz" # input fieldmap magnitude image - can be a 4D containing more than one
       FmapPhaseInputName="FieldMap_Phase.nii.gz" # input fieldmap phase image - in radians
       TE="2.46" #delta TE in ms for field map
       #Scan Settings for T1w
       T1wSampleSpacing="0.0000074" #DICOM field (0019,1018) in s
       UnwarpDir="z" 
       #Scan Settings for T2w
       T2wSampleSpacing="0.0000021" #DICOM field (0019,1018) in s

       ### If you use AvgrdcMethod="TOPUP" change variables in this part
       TopupConfig="NONE" #Topup config file

##############################
###### Functional input ######

#Naming conventions:	StudyFolder/SubjectFromSubjlist/task.nii.gz 
# A Tasklist is a list of names used to identify the particular fMRI experiments (these do not need to be the exact filenames but they need to be included in the filenames)
# Note that resting-state data is processed the same way as task data, so assign resting-state data to a "task" in the Tasklist

    # Write space delimited list of fMRI tasks in the Tasklist
    # No FMRI processing   
    Tasklist=""

    # The filename should use the wildcard %Tasklist% as part of the full filename
    # e.g. if the filename is "FUNC_r.nii.gz" (and "r" is in Tasklist) write the input as "FUNC_%TaskName%.nii.gz"
    fMRITimeSeries="%TaskName%.nii.gz"

    fMRISBRef="NONE" #A single band reference image (SBRef) is recommended if using multiband, set to NONE if you want to use the first volume of the timeseries for motion correction, or set it to the filename, including %TaskName% where appropriate
    FinalFMRIResolution="2" #Target final resolution of volumetric and surface fMRI data (in mm). 2mm is recommended.  
    SmoothingFWHM="2" #Smoothing on the surface; Recommended to be roughly the voxel size

   # The variables for distortion correction have "ForFUNC" in the name to distinguish them from the ones used for the structural scans
   DistortionCorrection="FIELDMAP" #FIELDMAP or TOPUP, distortion correction is required for accurate processing

       ### If you use DistortionCorrection="TOPUP" change variables in this part and read additional notes at the bottom of this file, otherwise skip it
        SpinEchoPhaseEncodeNegative="NONE" #For the spin echo field map volume with a negative phase encoding direction (LR in HCP data)
        SpinEchoPhaseEncodePositive="NONE" #For the spin echo field map volume with a positive phase encoding direction (RL in HCP data)
        TopUpConfigForFUNC="NONE" #Topup config if using TOPUP

       ### If you use DistortionCorrection="FIELDMAP" change variables in this part. Fieldmap can be used only for Siemens scanners.
         # If subject ID number is in the filename, replace it  with a wildcard %subjectID%
        FmapMagnitudeInputNameForFUNC="mag_raw_concat_%subjectID%.nii.gz" #Expects 4D Magnitude volume with two 3D timepoints
        FmapPhaseInputNameForFUNC="ph_MB_%subjectID%.nii.gz" #Expects a 3D Phase volume
        DwellTime="0.00072" # Effective Echo Spacing (or Dwelltime) of fMRI image in seconds - note that this is the echo spacing divided by any in-plane acceleration factor
        DeltaTE="2.46" #2.46ms for 3T, 1.02ms for 7T
        UnwarpdirForFUNC="-y"  
  
  
##########################################################################
######## Standard settings - do not normally need to be changed ##########


### Gradient non-linearity distortion correction (very important for HCP data!) ###
  # This can be left at NONE for most scanners except the HCP Skyra (or wherever the gradients are very non-linear or the isocentre is far from the brain centre)
  #GradientDistortionCoeffs="/home/fs0/stam/scratch/Pipelines/global/config/coeff_SC72C_Skyra.grad"
  GradientDistortionCoeffs="NONE"


### Miscellaneous configurations ###
  BrainSize="150" #BrainSize in mm, 150 for humans
  GrayordinatesResolution="2" #Usually 2mm
  HighResMesh="164" #Usually 164k vertices
  LowResMesh="32" #Usually 32k vertices
  

#### TEMPLATES AND LABELS LOCATIONS ####

  #Templates for T1w
  T1wTemplate="${HCPPIPEDIR_Templates}/MNI152_T1_0.7mm.nii.gz" #MNI0.7mm template
  T1wTemplateBrain="${HCPPIPEDIR_Templates}/MNI152_T1_0.7mm_brain.nii.gz" #Brain extracted MNI0.7mm template
  T1wTemplate2mm="${HCPPIPEDIR_Templates}/MNI152_T1_2mm.nii.gz" #MNI2mm template
  TemplateMask="${HCPPIPEDIR_Templates}/MNI152_T1_0.7mm_brain_mask.nii.gz" #Brain mask MNI0.7mm template
  Template2mmMask="${HCPPIPEDIR_Templates}/MNI152_T1_2mm_brain_mask_dil.nii.gz" #MNI2mm template

  FNIRTConfig="${HCPPIPEDIR_Config}/T1_2_MNI152_2mm.cnf" #FNIRT 2mm T1w Config

  # Templates for T2w (these can be ignored if no T2w images are being processed) 
  T2wTemplate="${HCPPIPEDIR_Templates}/MNI152_T2_0.7mm.nii.gz" #MNI0.7mm T2wTemplate
  T2wTemplateBrain="${HCPPIPEDIR_Templates}/MNI152_T2_0.7mm_brain.nii.gz" #Brain extracted MNI0.7mm T2wTemplate
  T2wTemplate2mm="${HCPPIPEDIR_Templates}/MNI152_T2_2mm.nii.gz" #MNI2mm T2wTemplate


  # Templates and labels for surface processing
  SurfaceAtlasDIR="${HCPPIPEDIR_Templates}/standard_mesh_atlases" #(Need to rename make surf.gii and add 32k)
  GrayordinatesSpaceDIR="${HCPPIPEDIR_Templates}/91282_Greyordinates" #(Need to copy these in)
  SubcorticalGrayLabels="${HCPPIPEDIR_Config}/FreeSurferSubcorticalLabelTableLut.txt"
  FreeSurferLabels="${HCPPIPEDIR_Config}/FreeSurferAllLut.txt"
  ReferenceMyelinMaps="${HCPPIPEDIR_Templates}/standard_mesh_atlases/Conte69.MyelinMap_BC.164k_fs_LR.dscalar.nii"



##############################
####### Detailed notes #######

#TOPUP

#To get accurate EPI distortion correction with TOPUP, the flags in PhaseEncodinglist must match the phase encoding
#direction of the EPI scan, and you must have used the correct images in SpinEchoPhaseEncodeNegative and Positive
#variables.  If the distortion is twice as bad as in the original images, flip either the order of the spin echo
#images or reverse the phase encoding list flag.  The pipeline expects you to have used the same phase encoding
#axis in the fMRI data as in the spin echo field map data (x/-x or y/-y).  
