import pandas as pd
import numpy as np

with open("assembly.fasta_rsem_index.idx.fa") as fp:
    lines = [i for i in fp.readlines()]

class RNAseq:
    def __init__(self) -> None:
        with open("gene_count_matrix.txt",mode =  "r",encoding="utf-8") as fp:
            header = [i.strip() for i in fp.readline().split("\t")]
            lines = [[float(e) if n != 0 and e != '' else e for n,e in enumerate(i.strip().split("\t"))] for i in fp.readlines()[:]]

        self.data = {
            'gene': [i[0] for i in lines],
            header[1] : [i[1] for i in lines],
            header[2] : [i[2] for i in lines],
            header[3] : [i[3] for i in lines],
            header[4] : [i[4] for i in lines],
            header[5] : [i[5] for i in lines],
            header[6] : [i[6] for i in lines],
            header[7] : [i[7] for i in lines],
            header[8] : [i[8] for i in lines],
            header[9] : [i[9] for i in lines],
            header[10] : [i[10] for i in lines],
            header[11] : [i[11] for i in lines],
            header[12] : [i[12] for i in lines],
            header[13] : [i[13] for i in lines],
            header[14] : [i[14] for i in lines],
            header[15] : [i[15] for i in lines],
            header[16] : [i[16] for i in lines],
            header[17] : [i[17] for i in lines],
            header[18] : [i[18] for i in lines],
        }

        with open("assembly.fasta_rsem_index.idx.fa") as fp:
            lines = [i for i in fp.readlines()] 
        genes: list[str] = []
        seqs:list[str] = []
        for n,i in enumerate(lines):
            if n%2 != 0:
                seqs.append(i.strip())
            else:
                genes.append(i.strip().split(" ")[0][1:20])
        self.bp_list: dict[str,str] = {i:j for i,j in zip(genes,seqs)
        }
    
    def analyse(self) -> None:
        df = pd.DataFrame(self.data)
        gene_lengths = {gene: len(seq) for gene, seq in self.bp_list.items()}
        df['gene_length'] = df['gene'].map(gene_lengths)

        # RPKMの計算（リード数 / 遺伝子長（kb））
        rpkm = df.iloc[:, 1:-1].div(df['gene_length'] / 1000, axis=0)

        # 各サンプルのRPKMの合計
        rpkm_sum = rpkm.sum()

        # TPMの計算（RPKM / RPKMの合計）* 1,000,000
        tpm = rpkm.div(rpkm_sum) * 1e6

        # 結果の出力
        print( tpm.head())  # 最初の数行を表示
        # set the fitst column as gene name
        tpm.index = df['gene']
        tpm.to_csv('tpm.csv', sep=',')  # tpmの出力

RNAseq().analyse()