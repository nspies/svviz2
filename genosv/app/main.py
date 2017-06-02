import logging
import sys
import time

from genosv.app import commandline
from genosv.app.datahub import DataHub
from genosv.visualize import visualize

logger = logging.getLogger(__name__)


def get_datahub():
    args = commandline.parse_args(sys.argv[1:])

    datahub = DataHub()
    datahub.set_args(args)
    datahub.align_distance = 0
    for sample_name, sample in datahub.samples.items():
        logger.info("Search distance: {:,}bp".format(sample.search_distance))
        datahub.align_distance = max(datahub.align_distance, sample.align_distance)
    logger.info("Align distance: {:,}bp".format(sample.align_distance))

    return datahub


def run(datahub):
    """ this runs the app on the provided datahub """
    for variant in datahub.get_variants():
        # if "86,868,786" not in str(variant):
            # print("skipping", str(variant))
            # continue
        t0 = time.time()
        temp_storage = datahub.genotype_cur_variant()
        t1 = time.time()
        print("TIME:::", t1-t0)
        visualize.visualize(datahub, temp_storage)

def main():
    """ entry point from command line """
    logging.basicConfig(level=logging.DEBUG)
    datahub = get_datahub()
    run(datahub)
