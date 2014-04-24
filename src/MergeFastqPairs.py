import sys
import itertools
import argparse

def main(parser):
    args = parser.parse_args()
    reads1 = args.reads1
    reads2 = args.reads2
    outfile = args.output
    print(args)

    i = 0
    for l12 in itertools.izip(reads1, reads2):
        l1, l2 = l12
        if i % 4 == 1 or i % 4 == 3:
            outfile.write(l1.rstrip())
            outfile.write(l2)
        else:
            outfile.write(l1)
        i += 1
        if i % 1000000 == 0:
            sys.stdout.write("\r\rprocessed {} reads".format(i))

    reads1.close()
    reads2.close()
    outfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge the paired read sequences in two files")
    parser.add_argument('reads1', nargs='?', type=argparse.FileType('r'))
    parser.add_argument('reads2', nargs='?', type=argparse.FileType('r'))
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'))
    main(parser)

