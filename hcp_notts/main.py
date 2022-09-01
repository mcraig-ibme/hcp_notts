#!/bin/env python
"""
HCP_NOTTS: Helper code for running HCP pipelines at Nottingham

Command line interface
"""
import os
import argparse
import sys
import logging

from ._version import __version__

LOG = logging.getLogger(__name__)

def makedirs(dpath, exist_ok=False):
    """
    Make directories, optionally ignoring them if they already exist
    """
    try:
        os.makedirs(dpath)
    except OSError as exc:
        import errno
        if not exist_ok or exc.errno != errno.EEXIST:
            raise

def setup_logging(logfile=None, log_level="info", log_stream=None):
    """
    Set the log level, formatters and output streams for the logging output

    By default this goes to <outdir>/logfile at level INFO
    """

    # Set log level on the root logger to allow for the possibility of 
    # debug logging on individual loggers
    level = getattr(logging, level.upper(), logging.INFO)
    logging.getLogger().setLevel(level)

    if logfile:
        # Send the log to an output logfile
        makedirs(os.path.dirname(logfile), True)
        logging.basicConfig(filename=logfile, filemode="w", level=level)

    if log_stream:
        # Can also supply a stream to send log output to as well (e.g. sys.stdout)
        extra_handler = logging.StreamHandler(log_stream)
        extra_handler.setFormatter(logging.Formatter('%(levelname)s : %(message)s'))
        logging.getLogger().addHandler(extra_handler)

def envdir(var):
    dpath = os.environ.get(var, "")
    if not dpath or not os.path.isdir(dpath):
        sys.stderr.write("{var} not set or is not a directory")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(f'HCP pipelines Nottingham helper script v{__version__}', add_help=True)
    g = parser.add_argument_group("Basic options")
    g.add_argument('--basedir', help='Path to directory containing subject subfolders', default=".")
    g.add_argument('--subjlist', help='Path to file containing list of subject subfolder names to run', required=True)
    g.add_argument('--t1w', help='Name of T1w file within subject subfolder. May contain {subjid} placeholder', default="T1w")
    g.add_argument('--t2w', help='Name of T2w file within subject subfolder, if available. May contain {subjid} placeholder')
    g.add_argument('--avgrdc-method', help='Averaging and readout distortion correction methods', default="NONE", choices=["NONE", "TOPUP", "FIELDMAP"])
    g.add_argument('--skip', nargs="*", help='Skip section of pipeline', choices=["PreFreeSurfer", "FreeSurfer", "PostFreeSurfer"])
    g.add_argument('--logfile', help='File to log script output')
    g.add_argument('--logdir', help='Name of subfolder for subject logs')

    g = parser.add_argument_group("Fieldmap options")
    g.add_argument('--fmap-mag', help='Fieldmap magnitude file name', default="FieldMap_Magnitude.nii.gz")
    g.add_argument('--fmap-phase', help='Fieldmap phase file name', default="FieldMap_Phase.nii.gz")
    g.add_argument('--fmap-te', help='Fieldmap TE', type=float, default=2.46)
    g.add_argument('--fmap-t1w-sample-spacing', help='Scan Settings for T1w: DICOM field (0019,1018) in s', type=float, default=0.0000074)
    g.add_argument('--fmap-t2w-sample-spacing', help='Scan Settings for T1w: DICOM field (0019,1018) in s', type=float, default=0.0000021)
    g.add_argument('--fmap-unwarp-dir', help='Unwarp direction', default='z')

    g = parser.add_argument_group("TOPUP options")
    g.add_argument('--topup-config', help='TOPUP config file')

    g = parser.add_argument_group("GDC options")
    g.add_argument('--gdc-coeffs', help='', default="NONE")

    g = parser.add_argument_group("Extended options")
    g.add_argument('--brain-size', help='BrainSize in mm, 150 for humans', type=float, default=150)
    g.add_argument('--grayordinates-res', help='Grayordinates resolution in mm', type=float, default=2)
    g.add_argument('--hires-mesh', help='High-resolution mesh (1000s of vertices)', type=int, default=164)
    g.add_argument('--lores-mesh', help='Low-resolution mesh (1000s of vertices)', type=int, default=32)

    g = parser.add_argument_group("Templates - non-absolute paths are relative to HCP template dir")
    g.add_argument('--template-t1w', help='T1w template', default="MNI152_T1_0.7mm.nii.gz")
    g.add_argument('--template-t1w-template-brain', help='T1w brain template', default="MNI152_T1_0.7mm_brain.nii.gz")
    g.add_argument('--template-mask', help='Template brain mask', default="MNI152_T1_0.7mm_brain_mask.nii.gz")
    g.add_argument('--template-t1w-2mm', help='T1w 2mm template', default="MNI152_T1_2mm.nii.gz")
    g.add_argument('--template-2mm-mask', help='2mm template brain mask', default="MNI152_T1_2mm_brain_mask_dil.nii.gz")
    g.add_argument('--template-t2w', help='T2w template', default="MNI152_T2_0.7mm.nii.gz")
    g.add_argument('--template-t2w-brain', help='T2w brain template', default="MNI152_T2_0.7mm_brain.nii.gz")
    g.add_argument('--template-t2w-2mm', help='T2w 2mm template', default="MNI152_T2_2mm.nii.gz")

    g = parser.add_argument_group("Config files - non-absolute paths are relative to HCP config dir")
    g.add_argument('--config-fnirt', help='', default="T1_2_MNI152_2mm.cnf")
    g.add_argument('--config-freesurfer-labels', help='', default="FreeSurferAllLut.txt")
    g.add_argument('--config-subcort-gray-labels', help='', default="FreeSurferSubcorticalLabelTableLut.txt")

    # FIXME below should be relative to templates dir???
    g.add_argument('--template-surface-atlas-dir', help='', default="standard_mesh_atlases")
    g.add_argument('--template-grayordinates-space-dir', help='', default="91282_Greyordinates")
    g.add_argument('--template-ref-myelin-maps', help='', default="standard_mesh_atlases/Conte69.MyelinMap_BC.164k_fs_LR.dscalar.nii")

    args = parser.parse_args()

    # Basic sanity checks
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(args.basedir):
        parser.error(f"Subject subfolder directory {args.basedir} does not exist")
    if not os.path.exists(args.subjlist):
        parser.error(f"Subject list file {args.subjlist} does not exist")

    setup_logging(args.logfile, "debug" if args.debug else "info", log_stream=sys.stdout)

    # For string formatting convenient to have arguments as dictionary
    argdict = vars(args)

    # Start up pipeline and get key 
    LOG.info(f"HCP Pipelines Nottingham wrapper v{__version__}")
    LOG.info("Subject folders located in: {basedir}".format(argdict))
    LOG.info("Subject list file: {subjlist}".format(argdict))

    argdict["fsldir"] = envdir("FSLDIR")
    argdict["hcpdir"] = envdir("HCPPIPEDIR")
    argdict["hcp_templ_dir"] = envdir("HCPPIPEDIR_Templates")
    argdict["hcp_config_dir"] = envdir("HCPPIPEDIR_Config")

    # If SGE_ROOT is set, use long.q FIXME how does this relate to Nottingham setup?
    argdict["queue"] = ""
    if "SGE_ROOT" in os.environ:
        argdict["queue"] = "-q long.q"

    # Apply standard prefix directories to templates and config dirs
    for k in argdict:
        if k.startswith("template") and not os.path.isabs(argdict[k]):
            argdict[k] = os.path.join(argdict["hcp_templ_dir"], argdict[k])
        if k.startswith("config") and not os.path.isabs(argdict[k]):
            argdict[k] = os.path.join(argdict["hcp_config_dir"], argdict[k])

    # Iterate through subjects, starting the pipeline jobs required
    with open(args.subjlist, "r") as f:
        for subjid in f:
            subjid = subjid.strip()
            subjdir = os.path.join(args.basedir, subjid)
            if not os.path.exists(subjdir):
                LOG.warn("Subject directory not found for {subjid} - skipping")
                continue

            t1fname = args.t1w.format(subjid=subjid)
            if not os.path.exists(t1fname):
                LOG.warn("T1w image: {t1fname} not found for {subjid} - skipping")
                continue

            if args.t2w:
                t2fname = args.t2w.format(subjid=subjid)
                if not os.path.exists(t2fname):
                    LOG.warn("T2w image: {t2fname} not found for {subjid} - skipping")
                    continue
            else:
                t2fname = "NONE"

            argdict["subjid"] = subjid
            argdict["subjid"] = subjdir
            argdict["subjlogdir"] = os.path.join(subjdir, args.logdir)
            argdict["t1fname"] = t1fname
            argdict["t2fname"] = t1fname
            makedirs(argdict["subjlogdir"], exist_ok=True)

            if "PreFreeSurfer" not in args.skip:
                # Run pre-freesurfer
                cmd = '{fsldir}/bin/fsl_sub {queue} \
                    -l {subjlogdir} \
                    {hcpdir}/PreFreeSurfer/PreFreeSurferPipeline.sh \
                    --path="{basedir}" \
                    --subject="{subjid}" \
                    --t1="{t1fname}" \
                    --t2="{t2fname}" \
                    --t1template="{template_t1w}" \
                    --t1templatebrain="{template_t1w_brain}" \
                    --t1template2mm="{template_t1w_2mm}" \
                    --t2template="{template_t2w}" \
                    --t2templatebrain="{template_t2w_brain}" \
                    --t2template2mm="{template_t2w_2mm}" \
                    --templatemask="{template_mask}" \
                    --template2mmmask="{template_mask_2mm}" \
                    --brainsize="{brain_size}" \
                    --fnirtconfig="{config_fnirt}" \
                    --avgrdcmethod="{avgrdc_method}" \
                    --fmapmag="{fmap_mag}" \
                    --fmapphase="{fmap_phase}" \
                    --echospacing="{fmap_te}" \
                    --t1samplespacing="{fmap_t1w_sample_spacing}" \
                    --t2samplespacing="{fmap_t2w_sample_spacing}" \
                    --unwarpdir="{fmap_unwarp_dir}" \
                    --gdcoeffs="{gdc_coeffs}" \
                    --topupconfig="{topup_config}"'.format(argdict)
                LOG.info(cmd)

            if "FreeSurfer" not in args.skip:
                # Run freesurfer
                # FIXME what if t2 does not exist?
                cmd = '{fsldir}/bin/fsl_sub {queue} \
                    -l {subjlogdir} \
                    -j $jobidPreFreeSurfer \
                    {hcpdir}/FreeSurfer/FreeSurferPipeline.sh \
                    --subject="{subjid}" \
                    --subjectDIR="{subjdir}/T1w" \
                    --t1="{subjdir}/T1w/T1w_acpc_dc_restore.nii.gz" \
                    --t1brain="{subjdir}/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
                    --t2="{subjdir}/T1w/T2w_acpc_dc_restore.nii.gz"'.format(argdict)
                LOG.info(cmd)
      
            if "PostFreeSurfer" not in args.skip:
                # Run post freesurfer
                # FIXME what if t2 does not exist?
                cmd = 'j{fsldir}/bin/fsl_sub {queue} \
                    -l {subjlogdir} \
                    -j $jobidFreeSurfer \
                    {hcpdir}/PostFreeSurfer/PostFreeSurferPipeline.sh \
                    --path="{basedir}" \
                    --subject="{subjid}" \
                    --surfatlasdir="{template_surface_atlas_dir}" \
                    --grayordinatesdir="{template_grayordinates_space_dir}" \
                    --grayordinatesres="{grayordinates_res}" \
                    --hiresmesh="{hires_mesh}" \
                    --lowresmesh="{lores_mesh}" \
                    --subcortgraylabels="{config_subcort_gray_labels}" \
                    --freesurferlabels="{config_freesurfer_labels}" \
                    --refmyelinmaps="{template_ref_myelin_maps}"'.format(argdict)
                LOG.info(cmd)
      
if __name__ == "__main__":
    main()
