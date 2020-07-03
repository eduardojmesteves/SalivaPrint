from argparse import ArgumentParser
from scipy.signal import find_peaks

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def build_args_parser():
    usage = 'python salivaanalysis.py\n       ' \
            'run with --help for arguments descriptions'
    parser = ArgumentParser(description='A python application that makes a statistical analysis over the SalivaTec'
                                        'database in order to use the extracted information to predict patients oral'
                                        'health.', usage=usage)
    parser.add_argument('-p', '--proteins', dest='proteins_dataset_path', default='data/proteins_dataset.csv',
                        help='Path to the proteins concentration dataset file. Should be in the CSV format.')
    parser.add_argument('-a', '--anamnesis', dest='anamnesis_dataset_path', default='data/anamnesis_dataset.csv',
                        help='Path to the anamnesis dataset file. Should be in the CSV format.')
    parser.add_argument('-o', '--out', dest='out_path', default='results',
                        help='Path to the directory where all the generated data will be saved.')

    return parser


def main():
    args_parser = build_args_parser()
    args = args_parser.parse_args()

    proteins_df = pd.read_csv(args.proteins_dataset_path)
    anamnesis_df = pd.read_csv(args.anamnesis_dataset_path)
    out_path = args.out_path

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    if not os.path.exists(out_path + '/peaks'):
        os.makedirs(out_path + '/peaks')

    # Normalizando os valores das concentrações proteicas em torno de 0
    proteins_df_norm = proteins_df.loc[:, proteins_df.columns != 'sample']
    proteins_size = len(proteins_df_norm.index)

    for i in range(0, proteins_size):
        min_value = proteins_df_norm.iloc[i].values.min()
        proteins_df_norm.iloc[i] -= min_value

    # Salvando informações dos picos
    x = proteins_df_norm.columns.values
    peaks_data = []

    for i in range(0, proteins_size):
        line = proteins_df_norm.iloc[i]

        plt.figure(figsize=(12, 6))

        plt.plot(x, line)

        plt.ylabel('Fluorescence Signal')
        plt.xlabel('Protein Molecular Weight (kDa)')
        plt.ylim(bottom=-20)
        plt.axhline(color='r')
        ticks, _ = plt.xticks()
        plt.xticks(ticks[::20])
        ticks, _ = plt.xticks()
        plt.xticks(ticks=ticks, labels=[j for j in x[::20]])

        peaks, properties = find_peaks(line, prominence=10, width=5, height=25)
        np.diff(peaks)

        plt.plot(peaks, line[peaks], "x")

        plt.vlines(x=peaks, ymin=line[peaks] - properties["prominences"], ymax=line[peaks], color="C1")
        plt.hlines(y=properties["width_heights"], xmin=properties["left_ips"], xmax=properties["right_ips"], color="C1")
        plt.savefig('results/peaks/' + proteins_df['sample'][i] + '_peaks.png')
        plt.close()

        # Calculando áreas sobre os picos
        cuts = (peaks[1:] + peaks[:-1]) // 2
        cuts = np.insert(cuts, [0, cuts.size], [0, x.size])
        peak_begins = np.zeros_like(peaks)
        peak_ends = np.zeros_like(peaks)
        areas = np.zeros(peaks.shape)
        peaks_values = np.zeros(peaks.shape)

        for j in range(peaks.size):
            peak_value = line.values[peaks[j]]
            y_cut = line.values[cuts[j]:cuts[j + 1]]
            baseline = np.median(y_cut)
            large = np.where(y_cut > 0.5 * (peak_value + baseline))[0]
            peak_begins[j] = large.min() + cuts[j]
            peak_ends[j] = large.max() + cuts[j]
            peaks_values[j] = peak_value
            areas[j] = np.sum(line.values[peak_begins[j]:peak_ends[j]] - baseline)

        peaks_data.append([proteins_df['sample'][i], peaks_values, properties['peak_heights'],
                           properties['prominences'], properties['width_heights'], areas])

    # Adicionando diagnóstico a planilha de proteínas
    proteins_df_norm.insert(0, 'sample', proteins_df['sample'])
    proteins_df_norm = pd.merge(proteins_df_norm, anamnesis_df[['sample', 'diagnosis']], on='sample', how='left')

    # Criando a planilha com os resultados dos picos
    peaks_df = pd.DataFrame(peaks_data, columns=['sample', 'peaks', 'heights', 'prominences', 'widths', 'areas'])
    peaks_df['peaks'] = [','.join(map(str, i)) for i in peaks_df['peaks']]
    peaks_df['heights'] = [','.join(map(str, i)) for i in peaks_df['heights']]
    peaks_df['prominences'] = [','.join(map(str, i)) for i in peaks_df['prominences']]
    peaks_df['widths'] = [','.join(map(str, i)) for i in peaks_df['widths']]
    peaks_df['areas'] = [','.join(map(str, i)) for i in peaks_df['areas']]
    peaks_df = pd.merge(peaks_df, anamnesis_df[['sample', 'diagnosis']], on='sample', how='left')

    # Salvando os resultados
    proteins_df_norm.to_csv('results/proteins_dataset_normalized.csv', index=False)
    peaks_df.to_csv('results/proteins_peaks_dataset.csv', index=False)


if __name__ == '__main__':
    main()