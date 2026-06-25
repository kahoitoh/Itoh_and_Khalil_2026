#!/usr/bin/env python3
import argparse
import pysam

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--in-bam", required=True)
    ap.add_argument("-m", "--map-tsv", required=True, help="readID<TAB>barcode")
    ap.add_argument("-o", "--out-bam", required=True)
    args = ap.parse_args()

    qname_to_cb = {}
    with open(args.map_tsv) as f:
        for line in f:
            qname, cb = line.rstrip("\n").split("\t")
            qname = qname.lstrip("@")
            qname_to_cb[qname] = cb

    n_total = n_tagged = n_unmatched = 0

    with pysam.AlignmentFile(args.in_bam, "rb") as fin, \
         pysam.AlignmentFile(args.out_bam, "wb", template=fin) as fout:
        for r in fin:
            n_total += 1
            cb = qname_to_cb.get(r.query_name)
            if cb is None:
                n_unmatched += 1
            else:
                r.set_tag("CB", cb, value_type="Z", replace=True)
                n_tagged += 1
            fout.write(r)

    print(f"total={n_total} tagged={n_tagged} unmatched={n_unmatched}")

if __name__ == "__main__":
    main()
