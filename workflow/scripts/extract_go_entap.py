import argparse
import re
import os

def main():
    
    parser = argparse.ArgumentParser(description='Extract GO terms from EnTap results')
    parser.add_argument('--entap_file', type=str,
                    help='file result from Entap')
    parser.add_argument('--kog_file', type=str,
                    help='file result from Entap')
    parser.add_argument('--type', type=str,
                    help='Type G for gene ou T for transcript')
    parser.add_argument('--output_go_per_gene', type=str,
                    help='output file')
    parser.add_argument('--output_gene_per_go', type=str,
                    help='output file')
    parser.add_argument('--output_species', type=str,
                    help='output file')
    parser.add_argument('--output_evalue', type=str,
                    help='output file')
    parser.add_argument('--output_kog', type=str,
                    help='output file')

    args = parser.parse_args()

    #print(args.entap_file, args.type, args.output)

    tlen = 5
    if args.type == 'G':
        tlen = 4
    
    go_out = dict()
    kog_out = dict()
    kog_names = dict()
    gene_out = dict()
    specie = dict()
    evalue_dist = {'1e-15 - 1e-5':0,'1e-30 - 1e-15':0,'1e-45 - 1e-30':0,'1e-60 - 1e-45':0,'1e-100 - 1e-60':0,'0 - 1e-100':0,'0':0}
    num_genes = 0
    func_list = [chr(x) for x in range(ord('A'), ord('Z') + 1)]
    func_list.remove('X')

    with open(args.kog_file, 'r') as kog:
        for line in kog:
            line=line.rstrip('\n\r')
            (kog_id, func) = line.split('\t')
            kog_names[kog_id]=func
    l=0
    with open(args.entap_file, 'r') as entap_file:
        entap_file.readline()
        for line in entap_file:
            line = line.rstrip('\n\r')
            locus_name = '_'.join(line.split('\t')[0].split('_')[0:tlen])

            specie[locus_name] = line.split('\t')[13] if line.split('\t')[13] != '' else 'None'
                
            evalue = line.split('\t')[10]
            if evalue != '':
                evalue = float(evalue)
                if evalue <= 1e-5 and evalue > 1e-15:
                    evalue_dist['1e-15 - 1e-5']+=1
                elif evalue <= 1e-15 and evalue > 1e-30:
                    evalue_dist['1e-30 - 1e-15']+=1
                elif evalue <= 1e-30 and evalue > 1e-45:
                    evalue_dist['1e-45 - 1e-30']+=1
                elif evalue <= 1e-45 and evalue > 1e-60:
                    evalue_dist['1e-60 - 1e-45']+=1
                elif evalue <= 1e-60 and evalue > 1e-100:
                    evalue_dist['1e-100 - 1e-60']+=1
                elif evalue <= 1e-100 and evalue > 0:
                    evalue_dist['0 - 1e-100']+=1
                elif evalue == 0:
                    evalue_dist['0']+=1

            kog = re.findall(r'KOG\d{4}',line.split('\t')[30])

            for kog_id in kog:
                if not locus_name in kog_out:
                    kog_out[locus_name]=set()
                kog_out[locus_name].add(kog_id)

            gos_bp = re.findall(r'GO:\d{7}-[\w|\s]+\(L=\d\),',line.split('\t')[18])
            gos_cc = re.findall(r'GO:\d{7}-[\w|\s]+\(L=\d\),',line.split('\t')[19])
            gos_mf = re.findall(r'GO:\d{7}-[\w|\s]+\(L=\d\),',line.split('\t')[20])

            gos_bp += re.findall(r'GO:\d{7}-[\w|\s]+\(L=\d\),',line.split('\t')[32])
            gos_cc += re.findall(r'GO:\d{7}-[\w|\s]+\(L=\d\),',line.split('\t')[33])
            gos_mf += re.findall(r'GO:\d{7}-[\w|\s]+\(L=\d\),',line.split('\t')[34])
                
            gos_bp.append('Biological Process')
            gos_cc.append('Celular Component')
            gos_mf.append('Molecular Function')

            gos = [gos_bp,gos_cc,gos_mf]
                
            #print(gos)
            #continue
            if not locus_name in go_out:
                go_out[locus_name]=set()
            for list_gos in gos:
                if len(list_gos) < 2:
                    continue
                go_type = list_gos.pop()
                for go in list_gos:
                    #g = gos = re.findall(r'GO:\d{7}',go)
                    #print(go)
                    #continue
                    (go_id, go_desc, level, _) = re.split('[-|(|)]', go)
                    #print(go_type, go_id, go_desc, level)
                        
                    if go_id in go_out[locus_name]:
                        continue
                    if not go_id in gene_out:
                        gene_out[go_id] = [go_type, go_desc, level,0]
                    gene_out[go_id][3]+=1
                    go_out[locus_name].add(go_id)
                        

    num_genes = len(go_out.keys())
    print(num_genes, l)
    with open(args.output_go_per_gene,'w') as fout:
        for g in go_out.keys():
            if len(go_out[g]) >= 1:
                fout.write(f"{g}\t{','.join(go_out[g])}\n")

    with open(args.output_gene_per_go,'w') as fout:
        for g in gene_out.keys():
            fout.write(f"{g}\t{gene_out[g][0]}\t{gene_out[g][3]}\t{gene_out[g][3]/num_genes}\t{gene_out[g][1]}\n")

    with open(args.output_species,'w') as fout:
        species = dict()
        for g in specie.keys():
            if not specie[g] in species:
                species[specie[g]]=0
            species[specie[g]]+=1

        for e in species.keys():
            fout.write(f"{e}\t{species[e]}\t{species[e]/num_genes}\n")

    with open(args.output_evalue,'w') as fout:
        for k in evalue_dist.keys():
            fout.write(f"{k}\t{evalue_dist[k]}\n")

    kog=dict()
    for l in func_list:
        kog[l]=0

    for gene in kog_out.keys():
        kogs=kog_out[gene]
        tmp=list()
        for k in kogs:
            if k in kog_names:
                if kog_names[k] in tmp or not kog_names[k] in func_list:
                    continue
                if not kog_names[k] in kog:
                    kog[kog_names[k]]=0
                kog[kog_names[k]]+=1
                tmp.append(kog_names[k])

    with open(args.output_kog, 'w') as fout:
        for k in sorted(kog.keys()):
            fout.write(f'{k}\t{kog[k]}\n')


    os.system(f'Rscript wego.R {args.output_gene_per_go}')
    os.system(f'Rscript pie_plot_specie.R {args.output_species}')
    os.system(f'Rscript pie_plot_evalue.R {args.output_evalue}')
    os.system(f'Rscript plot_cog.R {args.output_kog}')

    

if __name__ == "__main__":
    main()
