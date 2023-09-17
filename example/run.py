from argparse import ArgumentParser

from pyhfp import segment

from utils import load_vf, segments_split, visualize


def main(argv=None):
    """Parse arguments and run the given experiments."""

    parser = ArgumentParser(description="Tune executor")
    parser.add_argument(
        "-f",
        "--file-name",
        type=str,
        default="Mug.obj",
        help="File name of example to be segmented.",
    )
    parser.add_argument(
        "-n",
        "--num",
        type=int,
        default=7,
        help="Number of segments.",
    )

    args = parser.parse_args(argv)

    vs, fs = load_vf(args.file_name)
    idx = segment(vs, fs, args.num)
    segments = segments_split(vs, fs, idx)
    visualize(segments.values())
    

if __name__ == '__main__':
    main()
