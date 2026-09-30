# Copyright (C) 2013  Ian Harry
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
This module is responsible for setting up the matched-filtering stage of
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""


import os
import logging

from pycbc.workflow.core import FileList, make_analysis_dir
from pycbc.workflow.jobsetup import (select_matchedfilter_class,
                                     sngl_ifo_job_setup,
                                     multi_ifo_coherent_job_setup)

logger = logging.getLogger('pycbc.workflow.matched_filter')


def setup_matchedfltr_workflow(workflow, science_segs, datafind_outs,
                               tmplt_banks, output_dir=None,
                               injection_file=None, tags=None):
    '''
    This function aims to be the gateway for setting up a set of matched-filter
    jobs in a workflow. This function is intended to support multiple
    different ways/codes that could be used for doing this. For now the only
    supported sub-module is one that runs the matched-filtering by setting up
    a serious of matched-filtering jobs, from one executable, to create
    matched-filter triggers covering the full range of science times for which
    there is data and a template bank file.

    Parameters
    -----------
    Workflow : pycbc.workflow.core.Workflow
        The workflow instance that the coincidence jobs will be added to.
    science_segs : ifo-keyed dictionary of igwn_segments.segmentlist instances
        The list of times that are being analysed in this workflow.
    datafind_outs : pycbc.workflow.core.FileList
        An FileList of the datafind files that are needed to obtain the
        data used in the analysis.
    tmplt_banks : pycbc.workflow.core.FileList
        An FileList of the template bank files that will serve as input
        in this stage.
    output_dir : path
        The directory in which output will be stored.
    injection_file : pycbc.workflow.core.File, optional (default=None)
        If given the file containing the simulation file to be sent to these
        jobs on the command line. If not given no file will be sent.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names.

    Returns
    -------
    inspiral_outs : pycbc.workflow.core.FileList
        A list of output files written by this stage. This *will not* contain
        any intermediate products produced within this stage of the workflow.
        If you require access to any intermediate products produced at this
        stage you can call the various sub-functions directly.
    '''
    if tags is None:
        tags = []
    logger.info("Entering matched-filtering setup module.")
    make_analysis_dir(output_dir)
    cp = workflow.cp

    # Parse for options in .ini file
    mfltrMethod = cp.get_opt_tags("workflow-matchedfilter", "matchedfilter-method",
                                  tags)

    # Could have a number of choices here
    if mfltrMethod == "WORKFLOW_INDEPENDENT_IFOS":
        logger.info("Adding matched-filter jobs to workflow.")
        inspiral_outs = setup_matchedfltr_dax_generated(workflow, science_segs,
                                      datafind_outs, tmplt_banks, output_dir,
                                      injection_file=injection_file,
                                      tags=tags)
    elif mfltrMethod == "WORKFLOW_MULTIPLE_IFOS":
        logger.info("Adding matched-filter jobs to workflow.")
        inspiral_outs = setup_matchedfltr_dax_generated_multi(workflow,
                                      science_segs, datafind_outs, tmplt_banks,
                                      output_dir, injection_file=injection_file,
                                      tags=tags)
    else:
        errMsg = "Matched filter method not recognized. Must be one of "
        errMsg += "WORKFLOW_INDEPENDENT_IFOS or WORKFLOW_MULTIPLE_IFOS."
        raise ValueError(errMsg)

    logger.info("Leaving matched-filtering setup module.")
    return inspiral_outs

def setup_matchedfltr_dax_generated(workflow, science_segs, datafind_outs,
                                    tmplt_banks, output_dir,
                                    injection_file=None,
                                    tags=None):
    '''
    Setup matched-filter jobs that are generated as part of the workflow.
    This
    module can support any matched-filter code that is similar in principle to
    lalapps_inspiral, but for new codes some additions are needed to define
    Executable and Job sub-classes (see jobutils.py).

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    science_segs : ifo-keyed dictionary of igwn_segments.segmentlist instances
        The list of times that are being analysed in this workflow.
    datafind_outs : pycbc.workflow.core.FileList
        An FileList of the datafind files that are needed to obtain the
        data used in the analysis.
    tmplt_banks : pycbc.workflow.core.FileList
        An FileList of the template bank files that will serve as input
        in this stage.
    output_dir : path
        The directory in which output will be stored.
    injection_file : pycbc.workflow.core.File, optional (default=None)
        If given the file containing the simulation file to be sent to these
        jobs on the command line. If not given no file will be sent.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names.

    Returns
    -------
    inspiral_outs : pycbc.workflow.core.FileList
        A list of output files written by this stage. This *will not* contain
        any intermediate products produced within this stage of the workflow.
        If you require access to any intermediate products produced at this
        stage you can call the various sub-functions directly.
    '''
    if tags is None:
        tags = []
    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here ....

    cp = workflow.cp
    ifos = science_segs.keys()
    match_fltr_exe = os.path.basename(cp.get('executables','inspiral'))
    # Select the appropriate class
    exe_class = select_matchedfilter_class(match_fltr_exe)

    # Set up class for holding the banks
    inspiral_outs = FileList([])

    # Matched-filtering is done independently for different ifos, but might not be!
    # If we want to use multi-detector matched-filtering or something similar to this
    # it would probably require a new module
    for ifo in ifos:
        logger.info("Setting up matched-filtering for %s.", ifo)
        job_instance = exe_class(workflow.cp, 'inspiral', ifo=ifo,
                                               out_dir=output_dir,
                                               injection_file=injection_file,
                                               tags=tags)

        sngl_ifo_job_setup(workflow, ifo, inspiral_outs, job_instance,
                           science_segs[ifo], datafind_outs,
                           parents=tmplt_banks, allow_overlap=False)
    return inspiral_outs

def setup_matchedfltr_dax_generated_multi(workflow, science_segs, datafind_outs,
                                          tmplt_banks, output_dir,
                                          injection_file=None,
                                          tags=None):
    '''
    Setup matched-filter jobs that are generated as part of the workflow in
    which a single job reads in and generates triggers over multiple ifos.
    This module can support any matched-filter code that is similar in
    principle to pycbc_multi_inspiral, but for new codes some additions are
    needed to define Executable and Job sub-classes (see jobutils.py).

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    science_segs : ifo-keyed dictionary of igwn_segments.segmentlist instances
        The list of times that are being analysed in this workflow.
    datafind_outs : pycbc.workflow.core.FileList
        A FileList of the datafind files that are needed to obtain the
        data used in the analysis, and (if requested by the user) the vetoes
        File and (if requested by the user) the search sky-grid File.
    tmplt_banks : pycbc.workflow.core.FileList
        An FileList of the template bank files that will serve as input
        in this stage.
    output_dir : path
        The directory in which output will be stored.
    injection_file : pycbc.workflow.core.File, optional (default=None)
        If given the file containing the simulation file to be sent to these
        jobs on the command line. If not given no file will be sent.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names.

    Returns
    -------
    inspiral_outs : pycbc.workflow.core.FileList
        A list of output files written by this stage. This *will not* contain
        any intermediate products produced within this stage of the workflow.
        If you require access to any intermediate products produced at this
        stage you can call the various sub-functions directly.
    '''
    if tags is None:
        tags = []
    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here ....

    cp = workflow.cp
    ifos = sorted(science_segs.keys())
    match_fltr_exe = os.path.basename(cp.get('executables','inspiral'))

    # List for holding the output
    inspiral_outs = FileList([])

    logger.info("Setting up matched-filtering for %s.", ' '.join(ifos))

    if match_fltr_exe == 'pycbc_multi_inspiral':
        exe_class = select_matchedfilter_class(match_fltr_exe)
        bool_sg = ['make_sky_grid' in f.description for f in datafind_outs]
        n_sg = sum(bool_sg)
        if n_sg == 0:
            cp.set('inspiral', 'ra',
                   cp.get('workflow', 'ra'))
            cp.set('inspiral', 'dec',
                   cp.get('workflow', 'dec'))
        elif n_sg > 1:
            msg = f'{datafind_outs} has {n_sg} sky-grid files, '
            msg += 'instead of only one.'
            raise RuntimeError(msg)
        # Code lines for Fermi GBM are commented out for the time being
        # from pycbc.workflow.grb_utils import get_sky_grid_scale
        # if cp.has_option("jitter_skyloc", "apply-fermi-error"):
        #     cp.set('inspiral', 'sky-error',
        #            str(get_sky_grid_scale(float(cp.get('workflow',
        #                                                'sky-error')))))
        # else:
        #     cp.set('inspiral', 'sky-error',
        #            str(get_sky_grid_scale(float(cp.get('workflow',
        #                                                'sky-error')),
        #                                   sigma_sys=0.0)))
        # cp.set('inspiral', 'trigger-time',\
        #        cp.get('workflow', 'trigger-time'))
        # cp.set('inspiral', 'block-duration',
        #        str(abs(science_segs[ifos[0]][0]) - \
        #                2 * int(cp.get('inspiral', 'pad-data'))))

        job_instance = exe_class(workflow.cp, 'inspiral', ifo=ifos,
                                 out_dir=output_dir,
                                 injection_file=injection_file,
                                 tags=tags)
        if cp.has_option("workflow", "do-long-slides") and "slide" in tags[-1]:
            slide_num = int(tags[-1].replace("slide", ""))
            logger.info(
                "Setting up matched-filtering for slide %d",
                slide_num
            )
            slide_shift = int(cp.get("inspiral", "segment-length"))
            time_slide_dict = {ifo: (slide_num + 1) * ix * slide_shift
                               for ix, ifo in enumerate(ifos)}
            multi_ifo_coherent_job_setup(workflow, inspiral_outs, job_instance,
                                         science_segs, datafind_outs,
                                         output_dir, parents=tmplt_banks,
                                         slide_dict=time_slide_dict)
        else:
            multi_ifo_coherent_job_setup(workflow, inspiral_outs, job_instance,
                                         science_segs, datafind_outs,
                                         output_dir, parents=tmplt_banks)
    elif match_fltr_exe in ['pycbc_inspiral_fir', 'pycbc_inspiral']:
        exe_class = select_matchedfilter_class(match_fltr_exe)

        coincident_only = cp.has_option('workflow-matchedfilter', 'coincident-only')
        if coincident_only:
            analysis_segs = science_segs[ifos[0]]
            for ifo in ifos[1:]:
                analysis_segs = analysis_segs & science_segs[ifo]
        else:
            from igwn_segments import segmentlist
            analysis_segs = segmentlist([])
            for ifo in ifos:
                analysis_segs.extend(science_segs[ifo])
            analysis_segs.coalesce()

        if cp.has_option('workflow-matchedfilter', 'target-job-duration'):
            target_duration = int(cp.get('workflow-matchedfilter', 'target-job-duration'))
        else:
            target_duration = 2048

        if cp.has_option('workflow-matchedfilter', 'min-analysis-length'):
            min_analysis_length = int(cp.get('workflow-matchedfilter', 'min-analysis-length'))
        else:
            min_analysis_length = 0

        pad_data = 8
        if cp.has_option('inspiral', 'pad-data'):
            pad_data = int(cp.get('inspiral', 'pad-data'))

        from igwn_segments import segment

        job_instances = {}
        def get_job_instance(active_ifos):
            key = tuple(active_ifos)
            if key not in job_instances:
                job_instances[key] = exe_class(
                    workflow.cp, 'inspiral',
                    ifo=active_ifos if len(active_ifos) > 1 else active_ifos[0],
                    out_dir=output_dir,
                    injection_file=injection_file,
                    tags=tags
                )
            return job_instances[key]

        for curr_seg in analysis_segs:
            seg_len = abs(curr_seg)
            if min_analysis_length and seg_len < min_analysis_length:
                continue

            num_jobs = max(1, int(round(seg_len / float(target_duration))))
            job_dur = seg_len / float(num_jobs)

            for job_idx in range(num_jobs):
                v_start = int(round(curr_seg[0] + job_idx * job_dur))
                v_end = int(round(curr_seg[0] + (job_idx + 1) * job_dur))
                if job_idx == num_jobs - 1:
                    v_end = int(curr_seg[1])

                valid_seg = segment([v_start, v_end])
                data_seg = segment([v_start - pad_data, v_end + pad_data])

                if coincident_only:
                    job_ifos = list(ifos)
                else:
                    job_ifos = [ifo for ifo in ifos if science_segs[ifo].intersects_segment(valid_seg)]
                    if not job_ifos:
                        job_ifos = list(ifos)

                curr_job = get_job_instance(job_ifos)

                curr_dfouts = FileList([])
                if datafind_outs:
                    for ifo in job_ifos:
                        df = datafind_outs.find_all_output_in_range(ifo, data_seg, useSplitLists=True)
                        if df:
                            curr_dfouts.extend(df)

                for split_bank in tmplt_banks:
                    tag = list(tags)
                    if len(tmplt_banks) > 1 and hasattr(split_bank, 'tag_str') and split_bank.tag_str:
                        tag.append(split_bank.tag_str)
                    node = curr_job.create_node(
                        data_seg, valid_seg,
                        parent=split_bank,
                        df_parents=curr_dfouts if curr_dfouts else None,
                        tags=tag
                    )
                    workflow.add_node(node)
                    curr_out_files = [i for i in node.output_files if 'PSD_FILE' not in i.tags]
                    inspiral_outs.extend(curr_out_files)
    else:
        # Select the appropriate class
        raise ValueError("Not currently supported.")
    return inspiral_outs
