from lal import LIGOTimeGPS, YRJUL_SI
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from pycbc import version as pycbc_version
from pycbc import pnutils
from pycbc.tmpltbank import return_empty_sngl
import logging
import pycbc

class SingleCoincForGraceDB(object):
    """ Create xml files and submit them to gracedb from PyCBC Live """
    def __init__(self, ifos, coinc_results):
        # remember if this should be marked as HWINJ
        self.is_hardware_injection = False
        if 'foreground/HWINJ' in coinc_results:
            self.is_hardware_injection = True

        # Set up the bare structure of the xml document
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())

        proc_id = ligolw_process.register_to_xmldoc(outdoc, 'pycbc',
                     {}, ifos=ifos, comment='', version=pycbc_version.git_hash,
                     cvs_repository='pycbc/'+pycbc_version.git_branch,
                     cvs_entry_time=pycbc_version.date).process_id
    
        # Set up coinc_definer table
        coinc_def_table = lsctables.New(lsctables.CoincDefTable)
        coinc_def_id = lsctables.CoincDefID(0)
        coinc_def_row = lsctables.CoincDef()
        coinc_def_row.search = "inspiral"
        coinc_def_row.description = "sngl_inspiral<-->sngl_inspiral coincidences"
        coinc_def_row.coinc_def_id = coinc_def_id
        coinc_def_row.search_coinc_type = 0
        coinc_def_table.append(coinc_def_row)
        outdoc.childNodes[0].appendChild(coinc_def_table)

        # Set up coinc inspiral and coinc event tables
        coinc_id = lsctables.CoincID(0)
        coinc_event_table = lsctables.New(lsctables.CoincTable)
        coinc_event_row = lsctables.Coinc()
        coinc_event_row.coinc_def_id = coinc_def_id
        coinc_event_row.nevents = len(ifos)
        coinc_event_row.instruments = ','.join(ifos)
        coinc_event_row.time_slide_id = lsctables.TimeSlideID(0)
        coinc_event_row.process_id = proc_id
        coinc_event_row.coinc_event_id = coinc_id
        coinc_event_row.likelihood = 0.
        coinc_event_table.append(coinc_event_row)
        outdoc.childNodes[0].appendChild(coinc_event_table)

        # Get summary information for the coinc table
        mass1 = coinc_results['foreground/%s/mass1' % ifos[0]]
        mass2 = coinc_results['foreground/%s/mass1' % ifos[0]]
        mchirp, eta = pnutils.mass1_mass2_to_mchirp_eta(mass1, mass2)
        end_time = coinc_results['foreground/%s/end_time' % ifos[0]]
    
        # Set up the coinc inspiral table
        coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
        coinc_inspiral_row = lsctables.CoincInspiral()
        # This seems to be used as FAP, which should not be in gracedb
        coinc_inspiral_row.false_alarm_rate = 0
        coinc_inspiral_row.minimum_duration = 0.
        coinc_inspiral_row.set_ifos(ifos)
        coinc_inspiral_row.coinc_event_id = coinc_id
        coinc_inspiral_row.mchirp = mchirp
        coinc_inspiral_row.mass = mass1 + mass2
        coinc_inspiral_row.set_end(LIGOTimeGPS(end_time))
        coinc_inspiral_row.snr = coinc_results['foreground/stat']
        coinc_inspiral_row.combined_far = 1.0 / (YRJUL_SI * coinc_results['foreground/ifar'])
        coinc_inspiral_table.append(coinc_inspiral_row)
        outdoc.childNodes[0].appendChild(coinc_inspiral_table)

        # Set up sngls
        sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
        coinc_event_map_table = lsctables.New(lsctables.CoincMapTable)

        sngl_id = 0
        for ifo in ifos:
            names = [n.split('/')[-1] for n in coinc_results if 'foreground/%s' % ifo in n]
            sngl_id += 1
            sngl = return_empty_sngl()
            sngl.event_id = lsctables.SnglInspiralID(sngl_id)
            sngl.ifo = ifo
            for name in names:
                val = coinc_results['foreground/%s/%s' % (ifo, name)]
                if name == 'end_time':
                    sngl.set_end(LIGOTimeGPS(val))
                else:
                    try:
                        setattr(sngl, name, val)
                    except:
                        pass
            sngl.mtotal, sngl.eta = pnutils.mass1_mass2_to_mtotal_eta(
                    sngl.mass1, sngl.mass2)
            sngl.mchirp, junk = pnutils.mass1_mass2_to_mchirp_eta(
                    sngl.mass1, sngl.mass2)
            sngl.eff_distance = (sngl.sigmasq)**0.5 / sngl.snr
            sngl_inspiral_table.append(sngl)

            # Set up coinc_map entry
            coinc_map_row = lsctables.CoincMap()
            coinc_map_row.table_name = 'sngl_inspiral'
            coinc_map_row.coinc_event_id = coinc_id
            coinc_map_row.event_id = sngl.event_id
            coinc_event_map_table.append(coinc_map_row)

        outdoc.childNodes[0].appendChild(coinc_event_map_table)
        outdoc.childNodes[0].appendChild(sngl_inspiral_table)
        self.outdoc = outdoc

    def save(self, filename):
        ligolw_utils.write_filename(self.outdoc, filename)

    def upload(self, fname, psds, low_frequency_cutoff, testing=True):
#        from ligo.gracedb.rest import GraceDb
        import lal.series, lal
        if testing:
            group = 'Test'
        else:
            group = 'CBC'

        #r = gracedb.createEvent(group, "pycbc", fname, "AllSky").json()
        #logging.info("Uploaded event %s.", r["graceid"])    
        psds_lal = {}
        for ifo in psds:
            psd = psds[ifo]
            kmin = int(low_frequency_cutoff / psd.delta_f)
            fseries = lal.CreateREAL8FrequencySeries(
                "psd", psd.epoch, low_frequency_cutoff, psd.delta_f,
                lal.StrainUnit**2 / lal.HertzUnit, len(psd) - kmin)
            fseries.data.data = psd.numpy()[kmin:] / pycbc.DYN_RANGE_FAC ** 2.0
            psds_lal[ifo] = fseries
        
        psd_xmldoc = lal.series.make_psd_xmldoc(psds_lal)
        ligolw_utils.write_filename(psd_xmldoc, "tmp_psd.xml.gz", gz=True)

        #gracedb.writeLog(r["graceid"],
        #         "PyCBC PSD estimate from the time of event",
        #         "psd.xml.gz", open("tmp_psd.xml.gz", "rb").read(),
        #         "psd").json()
        #logging.info("Uploaded file psd.xml.gz to event %s.", r["graceid"])     



