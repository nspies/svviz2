import logging
import sys
import time

from genosv.app import commandline
from genosv.app.datahub import DataHub
from genosv.visualize import visualize
from genosv.app import report
from genosv.visualize import dotplots

FORMAT = '%(asctime)s - %(name)-25s - %(levelname)-5s - %(message)s'
DATEFMT = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
logger = logging.getLogger(__name__)


def get_datahub():
    args = commandline.parse_args(sys.argv[1:])

    datahub = DataHub()
    datahub.set_args(args)
    datahub.align_distance = 0
    for sample_name, sample in datahub.samples.items():
        logger.info("Search distance: {:,}bp".format(sample.search_distance))

    datahub.align_distance = max(sample.align_distance for sample in datahub.samples.values())
    if datahub.args.align_distance is not None:
        assert datahub.args.align_distance > 0, "--align-distance must be a positive integer"
        datahub.align_distance = datahub.args.align_distance
    logger.info("Align distance: {:,}bp".format(sample.align_distance))

    return datahub


def run(datahub):
    """ this runs the app on the provided datahub """
    for variant in datahub.get_variants():
        if not datahub.args.render_only:
            t0 = time.time()
            datahub.genotype_cur_variant()
            t1 = time.time()
            print("TIME:::", t1-t0)
            
        # visualize.visualize(datahub)
        # report.report(datahub)
        dotplots.generate_dotplots(datahub)
        
    datahub.cleanup()

def main():
    """ entry point from command line """
    logging.basicConfig(level=logging.DEBUG)
    datahub = get_datahub()
    run(datahub)
