import os
from rgt.GenomicRegionSet import GenomicRegionSet

"""
Creates a bed file with MPBSs with (in green) and without (in red) evidence.
Also, the name of the instances will be Y for evidence, N for non evidence.

Authors: Eduardo G. Gusmao, Zhijian Li
"""


def evidence_args(parser):
    # Input Options
    parser.add_argument("--mpbs-file", type=str, metavar="FILE", default=None,
                        help="motif predicted binding sites file. DEFAULT: None")
    parser.add_argument("--chip-file", type=str, metavar="FILE", default=None,
                        help="the ChIP-seq peak files. DEFAULT: None")
    parser.add_argument("--peak-ext", type=int, metavar="INT", default=100,
                        help="The number used to extend the ChIP-seq summit. DEFAULT: 100")

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written. DEFAULT: current directory")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default="evidence",
                        help="The prefix for results files. DEFAULT: evidence")


def evidence_run(args):
    # Expanding summits
    chip_summit_regions = GenomicRegionSet("TFBS Summit Regions")
    chip_summit_regions.read(args.chip_file)

    for region in iter(chip_summit_regions):
        summit = int(region.data.split()[-1]) + region.initial
        region.initial = max(summit - (args.peak_ext / 2), 0)
        region.final = summit + (args.peak_ext / 2)

    # Calculating intersections
    mpbs_regions = GenomicRegionSet("MPBS Regions")
    mpbs_regions.read(args.mpbs_file)

    chip_summit_regions.sort()
    mpbs_regions.sort()

    tfbs_regions = GenomicRegionSet("TFBS Regions")

    for mpbs_region in mpbs_regions:
        if chip_summit_regions.include(mpbs_region):
            mpbs_region.name = mpbs_region.name.split(":")[0] + ":Y"
        else:
            mpbs_region.name = mpbs_region.name.split(":")[0] + ":N"
        tfbs_regions.add(mpbs_region)

    tfbs_regions.sort()

    tfbs_fname = os.path.join(args.output_location, "{}.bed".format(args.output_prefix))
    tfbs_regions.write(tfbs_fname)
