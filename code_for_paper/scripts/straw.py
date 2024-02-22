import hicstraw
import argparse

def read_hic(hic_path):
    hic = hicstraw.HiCFile(hic_path)
    chrom_list= []

    for chrom in hic.getChromosomes():
        if (chrom.name != "All" and chrom.name != "MT" and chrom.name != "chrM"):
            chrom_list.append(chrom.name)

    chrom_list.reverse()
    chrom_list_str = ""

    for chrom in chrom_list:
        chrom_list_str += f"{chrom} "

    print(chrom_list_str)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='straw.py',
        allow_abbrev=False,
        description='What the program does', epilog='Text at the bottom of help'
    )
    parser.add_argument(
        "--hic_path",
        type=str,
        required=True,
        help="Input blablabla"
    )

    args = parser.parse_args()
    read_hic(hic_path=args.hic_path)