import pysam

with pysam.VariantFile("-", 'r') as VCF:

    # fix the header
    header = str(VCF.header).strip().splitlines()
    for line in header[:-1]:
        if "length=0" not in line:
            print(line)
    print("##INFO=<ID=IS_TRF,Number=.,Type=String,Description=\"\"")
    print("##INFO=<ID=SVCLASS,Number=.,Type=String,Description=\"\"")
    print("##INFO=<ID=CALLSET,Number=.,Type=String,Description=\"\"")
    print("##INFO=<ID=UNION,Number=.,Type=String,Description=\"\"")
    print(header[-1])

    # fix the "INFO/END" fields
    for variant in VCF:
        variant.stop = variant.pos - variant.info['SVLEN']
        print(str(variant).rstrip())

