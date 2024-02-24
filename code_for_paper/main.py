import argparse
from experiments import rao_2014
from experiments import lieberman_2009 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='main.py',
        allow_abbrev=False,
        description='What the program does',
        epilog='Text at the bottom of help'
    )
    parser.add_argument(
        "--docker_volume_path",
        type=str,
        required=True,
        help="docker_volume_path"
    )

    args = parser.parse_args()
    docker_volume_path = args.docker_volume_path

    rao_2014.run_all(docker_volume_path)
    # lieberman_2009.run_all(docker_volume_path)