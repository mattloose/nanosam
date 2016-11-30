import pysam
import h5py
from Bio import SeqIO
import string
import numpy as np
import sys,os
import argparse


class FileExist(argparse.Action):
    """Check if the input file exist."""
    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.exists(values):
             raise RuntimeError("File/path for '{}' does not exist, {}".format(self.dest, values))
        setattr(namespace, self.dest, values)


def get_parser():
    parser = argparse.ArgumentParser(
        description="""A simple tool to dump nanopore fast5s into a BAM type structure - a work in progress.\n Example command with supplied files:\n python nanosam_encode.py allreads.fastq allreads.sam test.bam\r\n This assumes the presence of fast5 files in the folder testpass. Command must be run from this folder. Running this command should generate a bam file equivalent to the ref.bam provided. """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("inputfastq", action=FileExist,
        help="Path to fastq file of reads extracted by poretools. nanosam assumes that the read files are in the same location as given in the poretools headers. If not, it will not be able to find the reads.")
    parser.add_argument("inputsam", action=FileExist,
        help="Path to a sam file with read alignment data. To do - check if this reads bam!")
    parser.add_argument("outputbam", type=str,
        help="Name of bam file to write out")
    #parser.add_argument("--threshold", default=90, type=int,
    #    help="Minimum match score to accept called barcodes.")

    return parser

def events_dict_to_numpy(d):
    events = np.empty(len(d["start"]), dtype=[('start', float), ('length', float),
                                            ('mean', float), ('stdv', float)])
    events["start"] = np.array(d["start"], dtype = float)
    events["length"] = np.array(d["length"], dtype = float)
    events["mean"] = np.array(d["mean"], dtype = float)
    events["stdv"] = np.array(d["stdv"], dtype = float)
    return events

def events_numpy_to_dict(events):
    d = {}
    try:
        d["start"] = events["start"].tolist()
        d["length"] = events["length"].tolist()
        d["mean"] = events["mean"].tolist()
        d["stdv"] = events["stdv"].tolist()
    except ValueError:
        d["start"] = events["start"].tolist()
        d["length"] = events["length"].tolist()
        d["mean"] = events["mean"].tolist()
        d["stdv"] = events["variance"].tolist()
    return d

def return_diff_list(listin):
    returnlist = list()
    for idx, val in enumerate(listin):
        if idx == 0:
            #print idx,val,val
            returnlist.append(val)
        else:
            #print idx,val,val-(listin[idx-1])
            returnlist.append(val-(listin[idx-1]))
    return returnlist

def main():

    args = get_parser().parse_args()

    samfile = pysam.AlignmentFile(args.inputsam, "r")
    header=dict()
    header["HD"]=dict()
    header["HD"]["SO"] = "unsorted"
    outsam = pysam.AlignmentFile(args.outputbam,"wb", template=samfile,header=header)
    fasta_sequences = SeqIO.parse(open(args.inputfastq),'fastq')

    fastalocationdict=dict()
    fastalendict=dict()

    for fasta in fasta_sequences:
            #print len(str(fasta.seq))
            name, sequence,description = fasta.id, str(fasta.seq),fasta.description
            #print name
            fastalocationdict[name]=description.split()[-1]
            fastalendict[name]=len(str(fasta.seq))
            #print description
            #print sequence

    counter = 0
    for samline in samfile.fetch():
        counter +=1
        print "processed read", counter, "\r",
        sys.stdout.flush()
        #print type(samline)
        #print samline.is_supplementary
        if samline.is_secondary is False and samline.is_supplementary is False:
            #print samline.infer_query_length(),fastalendict[samline.query_name]
            #print fastalocationdict[samline.query_name]
            hdf = h5py.File(fastalocationdict[samline.query_name], 'r')
            for read in hdf['Analyses']['EventDetection_001']['Reads']:
                events = hdf['Analyses']['EventDetection_001']['Reads'][read]['Events'][()]
            for read in hdf['Raw']['Reads']:
                signal = hdf['Raw']['Reads'][read]['Signal'][()]
            events = events_numpy_to_dict(events)

            listarray = signal.tolist()
            listarray = return_diff_list(listarray)
            startarray = return_diff_list(events["start"])
            #print startarray[0:100]
            lengtharray = return_diff_list(events["length"])
            #print lengtharray[0:100]
            #print type(events["mean"][0])
            samline.set_tag("oR", listarray)
            samline.set_tag("oS", startarray)
            samline.set_tag("oL", lengtharray)
            samline.set_tag("oM", events["mean"])
            samline.set_tag("oV", events["stdv"])
            hdf.close()
            outsam.write(samline)

    samfile.close()
    outsam.close()
    print "\nDone!"


if __name__ == "__main__":
    main()
