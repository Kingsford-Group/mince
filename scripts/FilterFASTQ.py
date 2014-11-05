import HTSeq
import sys

def main():
    ifile = sys.argv[1]
    ofile = sys.argv[2]
    r = HTSeq.FastqReader(ifile)
    with open(ofile, 'wb') as outstream:
        for ent in r:
            ent.write_to_fastq_file(outstream)

if __name__ == "__main__":
    main()
