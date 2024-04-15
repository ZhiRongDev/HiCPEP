import argparse
from experiments import rao_2014
from experiments import lieberman_2009 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='main.py',
        allow_abbrev=False,
        description='This program is used for producing the experiment results.',
        epilog='Text at the bottom of help'
    )
    parser.add_argument(
        "--data_store",
        type=str,
        required=True,
        help="This is the parameter specifying the path of the data_store directory."
    )

    args = parser.parse_args()
    data_store = args.data_store

    rao_2014.run_all(data_store)
    # lieberman_2009.run_all(data_store)