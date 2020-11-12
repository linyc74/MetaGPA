import pandas as pd
import ngslite as ngs
import matplotlib.pyplot as plt
from typing import Optional, List
from itertools import combinations


class GetUnalignedAminoacidPosition:

    faa: str
    header: str
    aligned_position: int

    def get_protein_seq(self) -> str:

        protein_seq = None

        with ngs.FastaParser(self.faa) as parser:
            for head, seq in parser:
                if head == self.header:
                    protein_seq = seq

        assert protein_seq is not None

        return protein_seq

    def get_unaligned_position(
            self, protein_seq: str) -> Optional[int]:

        pos = self.aligned_position - 1

        aa = protein_seq[pos]

        if aa == '-':
            return None

        aas_before = [a for a in protein_seq[:pos] if a != '-']

        return len(aas_before) + 1

    def main(
            self,
            faa: str,
            header: str,
            aligned_position: int) -> Optional[int]:
        """
        aligned_position: 1-based

        Returns 1-based position in the original unaligned protein sequence
        """

        self.faa = faa
        self.header = header
        self.aligned_position = aligned_position

        protein_seq = self.get_protein_seq()

        unaligned_position = self.get_unaligned_position(
            protein_seq=protein_seq)

        return unaligned_position


class DifferentialConservationScore:

    # link to BLOSUM80: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM80#L27
    BLOSUM = pd.read_csv('BLOSUM80/BLOSUM80.csv', index_col='id')

    aas_1: List[str]
    aas_2: List[str]

    def get_intra_group_similarity(self) -> float:

        scores = []

        for aas in [self.aas_1, self.aas_2]:

            c = 0  # count
            s = 0  # sum
            for a1, a2, in combinations(aas, 2):
                c += 1
                s += self.BLOSUM[a1][a2]

            scores.append(s / c)

        return sum(scores) / len(scores)

    def get_inter_group_similarity(self) -> float:
        c = 0  # count
        s = 0  # sum
        for a1 in self.aas_1:
            for a2 in self.aas_2:
                c += 1
                s += self.BLOSUM[a1][a2]
        return s / c

    def main(
            self,
            aminoacids1: List[str],
            aminoacids2: List[str]) -> float:

        self.aas_1 = [a for a in aminoacids1 if a != '-']
        self.aas_2 = [a for a in aminoacids2 if a != '-']

        a = self.get_intra_group_similarity()
        b = self.get_inter_group_similarity()

        return a - b


def fraction_present(aminoacids: List[str]) -> float:
    c = 0
    for a in aminoacids:
        if a != '-':
            c += 1
    return c / len(aminoacids)


class Main:

    aligned_faa = '../20_0904_thymidylate_synthase_muscle/20_0904_thymidylate_synthase_muscle.faa'

    category_1 = 'Modified'
    category_2 = 'Unmodified'

    header_prefix_to_category = {
        'E_coli':  'Unmodified',
        'NEBs1': 'Modified',
        'Modified': 'Modified',
        'Unmodified': 'Unmodified',
    }

    fraction_cutoff = 0.5

    ecoli_header = 'E_coli_Thymidylate_synthase'

    residue_data_csv = f'{__file__[:-3]}_residue_data.csv'
    output_csv = f'{__file__[:-3]}.csv'
    output_png = f'{__file__[:-3]}.png'

    def __get_category(self, header: str) -> Optional[str]:
        for key, val in self.header_prefix_to_category.items():
            if header.startswith(key):
                return val

    def get_residue_data_df(self) -> pd.DataFrame:

        data = {
            'category': [],
            'fasta_header': [],
            'position': [],
            'aminoacid': [],
        }

        with ngs.FastaParser(self.aligned_faa) as parser:
            for head, seq in parser:
                for pos, aa in enumerate(seq):

                    category = self.__get_category(header=head)

                    data['category'].append(category)
                    data['fasta_header'].append(head)
                    data['position'].append(pos+1)  # 0 to 1-based
                    data['aminoacid'].append(aa)

        return pd.DataFrame(data=data)

    def get_position_score_df(
            self,
            residue_data_df: pd.DataFrame) -> pd.DataFrame:

        data = {
            'position': [],
            'fraction_present': [],
            'differential_conservation_score': [],
        }

        df = residue_data_df
        positions = sorted(df['position'].unique())

        for pos in positions:

            all_aa = df.loc[df['position'] == pos, 'aminoacid']
            fraction = fraction_present(aminoacids=all_aa)
            if fraction < self.fraction_cutoff:
                continue

            aminoacids1 = df.aminoacid[(df['position'] == pos) & (df['category'] == self.category_1)]
            aminoacids2 = df.aminoacid[(df['position'] == pos) & (df['category'] == self.category_2)]

            score = DifferentialConservationScore().main(
                aminoacids1=aminoacids1,
                aminoacids2=aminoacids2)

            data['position'].append(pos)
            data['fraction_present'].append(fraction)
            data['differential_conservation_score'].append(score)

        return pd.DataFrame(data=data)

    def __to_ecoli_position(
            self, aligned_position: int) -> int:

        return GetUnalignedAminoacidPosition().main(
            faa=self.aligned_faa,
            header=self.ecoli_header,
            aligned_position=aligned_position)

    def add_ecoli_position(
            self, position_score_df: pd.DataFrame) -> pd.DataFrame:

        df = position_score_df.copy()

        df['ecoli_position'] = df['position'].apply(self.__to_ecoli_position)

        return df

    def plot_score(self, position_score_df: pd.DataFrame):

        df = position_score_df

        plt.figure(figsize=(6, 6))
        plt.plot(
            df['position'],
            df['differential_conservation_score'],
            marker='o',
            markersize=3,
            linestyle='None',
            markeredgewidth=0,
            markeredgecolor='black'
        )
        plt.xlabel('Position')
        plt.ylabel('Differential conservation score')
        plt.savefig(self.output_png, dpi=600)
        plt.close()

    def main(self):

        residue_data_df = self.get_residue_data_df()

        residue_data_df.to_csv(self.residue_data_csv, index=False)

        position_score_df = self.get_position_score_df(
            residue_data_df=residue_data_df)

        position_score_df = self.add_ecoli_position(
            position_score_df=position_score_df)

        position_score_df.to_csv(self.output_csv, index=False)

        self.plot_score(position_score_df=position_score_df)


if __name__ == '__main__':
    Main().main()
