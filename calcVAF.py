fileN = '/Users/valdezkm/Documents/LP/Val_ClinOmics_Compass_026_T1D_E.strelka.snvs.raw.vcf'

outFile = open('/Users/valdezkm/Documents/LP/gitHub/kv/vaf_py.vcf','w')

with open(fileN) as myfile:
    info = []
    notInfo = []
    chrom = []
    dat = []
    for each in myfile:                                                     # for each line in file
        if each.startswith('##'):                                           # subset commented section
            if each.startswith('##INFO'):                                   # subset info lines, will add vaf later
                info.append(each.strip())
            if not each.startswith('##INFO'):                               # if not info lines, subset
                notInfo.append(each.strip())
        if each.startswith('#CHROM'):                                       # pull out table header
            chrom.append(each.strip())
        if not each.startswith('#'):                                       # subset table and parse ACGTs, use tier 1
            As = float(each.split('\t')[10].split(':')[4].split(',')[0])
            Cs = float(each.split('\t')[10].split(':')[5].split(',')[0])
            Gs = float(each.split('\t')[10].split(':')[6].split(',')[0])
            Ts = float(each.split('\t')[10].split(':')[7].split(',')[0])

            REF = each.split('\t')[3]                                      # Find REF variant, set appropriate value for depth
            if REF == 'A':
                DEPTH = As
            elif REF == 'C':
                DEPTH = Cs
            elif REF == 'G':
                DEPTH = Gs
            elif REF == 'T':
                DEPTH = Ts

            ALT = each.split('\t')[4]                                       # find ALT variant, calc VAF
            if ALT == 'A':
                VAF = round(As / (As + DEPTH), 2)
            elif ALT == 'C':
                VAF = round(Cs / (Cs + DEPTH), 2)
            elif ALT == 'G':
                VAF = round(Gs / (Gs + DEPTH), 2)
            elif ALT == 'T':
                VAF = round(Ts / (Ts + DEPTH), 2)

            sepEach = each.split('\t')
            sepEach[7] = sepEach[7] + ';VAF=' + str(VAF)                    # new info column with vaf
            eachStr = '\t'
            each = eachStr.join(sepEach)                                    # redefine orig each variable
            dat.append(each.strip())

    info.append('##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">')   # add info tag vaf to header

    together = notInfo + info + chrom + dat                                                     # put vcf file back together
    for x in together:
        outFile.write(x.strip() + '\n')

