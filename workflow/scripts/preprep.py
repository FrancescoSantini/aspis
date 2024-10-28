import sys

def preprep(input_csv, input_gtf, output_lst, output_csv):
    with open(output_lst, 'w') as out0, open(output_csv, 'w') as out1:
        hdr = 'sample\tcondition\n'
        out1.write(hdr)
        with open(input_csv, 'r') as hnd:
            bioproject = input_gtf.split('/')[-2]  # Extract bioproject from input_gtf path
            for line in hnd.readlines()[1:]:
                smp, prj, cnd = line.strip().split(',')[1:4]
                if prj == bioproject:
                    oln0 = f"{smp}\tresults/assembly/{prj}/{smp}.gtf\n"
                    out0.write(oln0)
                    oln1 = f"{smp}\t{cnd}\n"
                    out1.write(oln1)

if __name__ == "__main__":
    input_csv, input_gtf, output_lst, output_csv = sys.argv[1:]
    preprep(input_csv, input_gtf, output_lst, output_csv)
