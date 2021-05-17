import argparse
import re
import os

def main():
    
    parser = argparse.ArgumentParser(description='Extract GO terms from EnTap results')
    parser.add_argument('--entap_file', type=str,
                    help='file result from Entap')
    parser.add_argument('--type', type=str,
                    help='Type G for gene ou T for transcript')
    parser.add_argument('--up', type=str,
                    help='Type G for gene ou T for transcript')
    parser.add_argument('--down', type=str,
                    help='Type G for gene ou T for transcript')

    args = parser.parse_args()

    tlen = 5
    if args.type == 'G':
        tlen = 4
        
    up_list = list()
    down_list = list()

    up_out = set()
    down_out = set()

    with open(args.up, 'r') as up:
        for line in up:
            line = line.rstrip('\n\r')
            up_list.append(line)

    with open(args.down, 'r') as down:
        for line in down:
            line = line.rstrip('\n\r')
            down_list.append(line)

    
    with open(args.entap_file, 'r') as entap_file:
        entap_file.readline()
        for line in entap_file:
            line = line.rstrip('\n\r')
            locus_name = '_'.join(line.split('\t')[0].split('_')[0:tlen])
            #if 'actinopterygii' in line:
            if True:
                if locus_name in up_list:
                    print(f"UP:{line}")
                    if line.split('\t')[27] != '':
                        up_out.add(line.split('\t')[27])
                elif locus_name in down_list:
                    if line.split('\t')[27] != '':
                        down_out.add(line.split('\t')[27])
    with open('up_enz_out2.txt', 'w') as fout:
        for gene in up_out:
            fout.write(f"{gene}\n")

    with open('down_enz_out2.txt', 'w') as fout:
        for gene in down_out:
            fout.write(f"{gene}\n")

if __name__ == "__main__":
    main()
