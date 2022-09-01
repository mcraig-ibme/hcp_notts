#!/bin/env python
"""
HCP_NOTTS: Helper code for running HCP pipelines at Nottingham

Command line interface
"""
import os
import argparse
import sys
import subprocess

from ._version import __version__

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

def envdir(var):
    dpath = os.environ.get(var, "")
    if not dpath or not os.path.isdir(dpath):
        sys.stderr.write(f"{var} not set or is not a directory\n")
        sys.exit(1)
    return dpath
  
def submit(cmd, argsdict, dep_job_id=None, queue=None):
    sub_cmd = ['{fsldir}/bin/fsl_sub', '-l', '{subjlogdir}']
    if queue:
        sub_cmd += ['-q', queue]
    if dep_job_id:
        sub_cmd += ['-j', str(dep_job_id)]
    sub_cmd += cmd
    sub_cmd = [s.format(**argsdict) for s in sub_cmd]
    if argsdict["debug"]:
      print(" ".join(sub_cmd))
    output = subprocess.check_output(sub_cmd)
    return int(output)

def main():
    parser = argparse.ArgumentParser(f'HCP pipelines Nottingham helper script v{__version__}', add_help=True)
    g = parser.add_argument_group("Basic options")
    g.add_argument('--basedir', help='Path to directory containing subject subfolders', default=".")
    g.add_argument('--subjlist', help='Path to file containing list of subject subfolder names to run', required=True)
    g.add_argument('--t1w', help='Name of T1w file within subject subfolder. May contain {subjid} placeholder', default="T1w")
    g.add_argument('--t2w', help='Name of T2w file within subject subfolder, if available. May contain {subjid} placeholder')
    g.add_argument('--avgrdc-method', help='Averaging and readout distortion correction methods', default="NONE", choices=["NONE", "TOPUP", "FIELDMAP"])
    g.add_argument('--skip', nargs="*", help='Skip section of pipeline', choices=["PreFreeSurfer", "FreeSurfer", "PostFreeSurfer"], default=[])
    g.add_argument('--logfile', help='File to log script output')
    g.add_argument('--logdir', help='Name of subfolder for subject logs', default='hcp_pipeline_logs')
    g.add_argument('--debug', help='Enable debug logging', action="store_true", default=False)

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
    g.add_argument('--template-t1w-brain', help='T1w brain template', default="MNI152_T1_0.7mm_brain.nii.gz")
    g.add_argument('--template-mask', help='Template brain mask', default="MNI152_T1_0.7mm_brain_mask.nii.gz")
    g.add_argument('--template-t1w-2mm', help='T1w 2mm template', default="MNI152_T1_2mm.nii.gz")
    g.add_argument('--template-mask-2mm', help='2mm template brain mask', default="MNI152_T1_2mm_brain_mask_dil.nii.gz")
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

    # For string formatting convenient to have arguments as dictionary
    argdict = vars(args)

    # Start up pipeline and get key 
    print(f"HCP Pipelines Nottingham wrapper v{__version__}")
    print("Subject folders located in: {basedir}".format(**argdict))
    print("Subject list file: {subjlist}".format(**argdict))

    argdict["fsldir"] = envdir("FSLDIR")
    argdict["hcpdir"] = envdir("HCPPIPEDIR")
    argdict["hcp_templ_dir"] = envdir("HCPPIPEDIR_Templates")
    argdict["hcp_config_dir"] = envdir("HCPPIPEDIR_Config")
    print(argdict)

    # If SGE_ROOT is set, use long.q FIXME how does this relate to Nottingham setup?
    if os.environ.get("SGE_ROOT", ""):
        queue = "long.q"
    else:
        queue = ""

    # Apply standard prefix directories to templates and config dirs
    for k in argdict:
        print(k)
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
                print(f"WARNING: Subject directory not found for {subjid} - skipping")
                continue

            t1fname = os.path.join(subjdir, args.t1w.format(subjid=subjid))
            if not os.path.exists(t1fname):
                print(f"WARNING: T1w image: {t1fname} not found for {subjid} - skipping")
                continue

            if args.t2w:
                t2fname = os.path.join(subjdir, args.t2w.format(subjid=subjid))
                if not os.path.exists(t2fname):
                    print(f"WARNING: T2w image: {t2fname} not found for {subjid} - skipping")
                    continue
            else:
                t2fname = "NONE"

            argdict["subjid"] = subjid
            argdict["subjdir"] = subjdir
            argdict["subjlogdir"] = os.path.join(subjdir, args.logdir)
            argdict["t1fname"] = t1fname
            argdict["t2fname"] = t1fname
            makedirs(argdict["subjlogdir"], exist_ok=True)

            if "PreFreeSurfer" not in args.skip:
                print("\nSubmitting pre-freesurfer")
                cmd = [
	            '{hcpdir}/PreFreeSurfer/PreFreeSurferPipeline.sh',
                    '--path="{basedir}"',
                    '--subject="{subjid}"',
                    '--t1="{t1fname}"',
                    '--t2="{t2fname}"',
                    '--t1template="{template_t1w}"',
                    '--t1templatebrain="{template_t1w_brain}"',
                    '--t1template2mm="{template_t1w_2mm}"',
                    '--t2template="{template_t2w}"',
                    '--t2templatebrain="{template_t2w_brain}"',
                    '--t2template2mm="{template_t2w_2mm}"',
                    '--templatemask="{template_mask}"',
                    '--template2mmmask="{template_mask_2mm}"',
                    '--brainsize="{brain_size}"',
                    '--fnirtconfig="{config_fnirt}"',
                    '--avgrdcmethod="{avgrdc_method}"',
                    '--fmapmag="{fmap_mag}"',
                    '--fmapphase="{fmap_phase}"',
                    '--echospacing="{fmap_te}"',
                    '--t1samplespacing="{fmap_t1w_sample_spacing}"',
                    '--t2samplespacing="{fmap_t2w_sample_spacing}"',
                    '--unwarpdir="{fmap_unwarp_dir}"',
                    '--gdcoeffs="{gdc_coeffs}"',
                    '--topupconfig="{topup_config}"',
                ]
                job_id_prefreesurfer = submit(cmd, argdict, queue=queue)

            if "FreeSurfer" not in args.skip:
                print("\nSubmitting freesurfer")
                # FIXME what if t2 does not exist?
                cmd = [
                    '{hcpdir}/FreeSurfer/FreeSurferPipeline.sh',
                    '--subject="{subjid}"',
                    '--subjectDIR="{subjdir}/T1w"',
                    '--t1="{subjdir}/T1w/T1w_acpc_dc_restore.nii.gz"',
                    '--t1brain="{subjdir}/T1w/T1w_acpc_dc_restore_brain.nii.gz"',
                    '--t2="{subjdir}/T1w/T2w_acpc_dc_restore.nii.gz"',
                ]
                job_id_freesurfer = submit(cmd, argdict, job_id_prefreesurfer, queue=queue)
      
            if "PostFreeSurfer" not in args.skip:
                print("\nSubmitting post-freesurfer")
                # FIXME what if t2 does not exist?
                cmd = [
                    '{hcpdir}/PostFreeSurfer/PostFreeSurferPipeline.sh',
                    '--path="{basedir}"',
                    '--subject="{subjid}"',
                    '--surfatlasdir="{template_surface_atlas_dir}"',
                    '--grayordinatesdir="{template_grayordinates_space_dir}"',
                    '--grayordinatesres="{grayordinates_res}"',
                    '--hiresmesh="{hires_mesh}"',
                    '--lowresmesh="{lores_mesh}"',
                    '--subcortgraylabels="{config_subcort_gray_labels}"',
                    '--freesurferlabels="{config_freesurfer_labels}"',
                    '--refmyelinmaps="{template_ref_myelin_maps}"',
                ]
                job_id_postfreesurfer = submit(cmd, argdict, job_id_freesurfer, queue=queue)

if __name__ == "__main__":
    main()
