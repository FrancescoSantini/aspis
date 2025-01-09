import sys

def preprep(input_csv, input_gtf, output_lst, output_csv):
    with open(output_lst, 'w') as out0, open(output_csv, 'w') as out1:
        # Update the header to include 'biosample'
        hdr = 'sample\tcondition\tbiosample\n'
        out1.write(hdr)
        with open(input_csv, 'r') as hnd:
            # Extract bioproject from input_gtf path
            bioproject = input_gtf.split('/')[-2]
            for line in hnd.readlines()[1:]:  # Skip the header row
                # Extract columns: sample_name, run, biosample, bioproject, condition
                smp, run, bsm, prj, cnd = line.strip().split(',')
                if prj == bioproject:
                    # Write to the GTF list file (unchanged)
                    oln0 = f"{run}\tresults/assembly/{prj}/{run}.gtf\n"
                    out0.write(oln0)
                    # Write to phenodata CSV, including biosample
                    oln1 = f"{run}\t{cnd}\t{bsm}\n"
                    out1.write(oln1)

if __name__ == "__main__":
    input_csv, input_gtf, output_lst, output_csv = sys.argv[1:]
    preprep(input_csv, input_gtf, output_lst, output_csv)
