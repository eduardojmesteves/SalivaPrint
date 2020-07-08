from argparse import ArgumentParser
from scipy.signal import find_peaks

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os


def build_args_parser():
    usage = 'python salivaanalysis.py\n       ' \
            'run with --help for arguments descriptions'
    parser = ArgumentParser(description='A python application that generates input data and makes a statistical'
                                        'analysis over the SalivaTec database in order to use the extracted information'
                                        'to predict patients oral health.', usage=usage)
    parser.add_argument('-p', '--proteins', dest='proteins_dataset_path', default='data/proteins_dataset.csv',
                        help='Path to the proteins concentration dataset file. Should be in the CSV format.')
    parser.add_argument('-a', '--anamnesis', dest='anamnesis_dataset_path', default='data/anamnesis_dataset.csv',
                        help='Path to the anamnesis dataset file. Should be in the CSV format.')
    parser.add_argument('-o', '--out', dest='out_path', default='results',
                        help='Path to the directory where all the generated data will be saved.')
    parser.add_argument('-c', '--analysis_choice', dest='analysis_choice', required=True,
                        choices=['gender', 'age_group', 'oral_diagnosis'], help='Type of analysis to be made')
    parser.add_argument('--save_graphics', action='store_true',
                        help='Use this flag to save graphics generated during the process')

    return parser


def do_levene(groups, analysis_variable_values, variable_name, analysis_variable_name):
    levene_results = []
    for a in analysis_variable_values:
        groups_aux = [group[a] for group in groups]
        levene_result = stats.levene(*groups_aux)
        levene_results.append(levene_result)

    levene_test = pd.DataFrame(levene_results)
    levene_test[analysis_variable_name] = analysis_variable_values
    levene_test = levene_test[[analysis_variable_name, 'statistic', 'pvalue']].round(3)
    levene_test.to_csv('results/' + variable_name + '_' + analysis_variable_name + '_levene_analysis.csv',
                       index=False)

    levene_test_result = levene_test.loc[levene_test['pvalue'] >= 0.05]

    if len(levene_test_result.index) >= 0.95 * len(levene_test.index):
        print('The levine test showed that the samples are homogeneous between ' + variable_name + ' variable with ' +
              analysis_variable_name)
        levene_test_result.to_csv('results/' + variable_name + '_mw_levene_analysis_95.csv', index=False)
    else:
        print('The levine test showed that the samples are NOT homogeneous between ' + variable_name + ' variable with '
              + analysis_variable_name)


def do_normality(group, analysis_variable_values, variable_name, group_name, analysis_variable_name):
    normality_results = []
    for a in analysis_variable_values:
        normality = stats.normaltest(group[a])
        normality_results.append(normality)

    normality_test = pd.DataFrame(normality_results)
    normality_test[analysis_variable_name] = analysis_variable_values
    normality_test = normality_test[[analysis_variable_name, 'statistic', 'pvalue']].round(3)

    normality_test.to_csv('results/' + variable_name + '_' + group_name + '_' + analysis_variable_name +
                          '_normality_analysis.csv', index=False)

    normality_test_result = normality_test.loc[normality_test['pvalue'] >= 0.05]

    if len(normality_test_result.index) >= 0.95 * len(normality_test.index):
        print('The normality test showed that the samples are homogeneous between ' + variable_name + ' "' + group_name
              + '" variable with ' + analysis_variable_name)
        normality_test_result.to_csv('results/' + variable_name + '_' + group_name + '_' + analysis_variable_name +
                                     '_normality_analysis_95.csv', index=False)
    else:
        print('The normality test showed that the samples are NOT homogeneous between ' + variable_name + ' "' +
              group_name + '" variable with ' + analysis_variable_name)


def do_mann_whitney(groups, analysis_variable_values, variable_name, analysis_variable_name):
    t_results = []
    for a in analysis_variable_values:
        groups_aux = [group[a] for group in groups]
        t_result = stats.mannwhitneyu(*groups_aux, use_continuity=True, alternative='two-sided')
        t_results.append(t_result)

    t_test = pd.DataFrame(t_results)
    t_test[analysis_variable_name] = analysis_variable_values
    t_test = t_test[[analysis_variable_name, 'statistic', 'pvalue']].round(3)
    t_test.to_csv('results/' + variable_name + '_' + analysis_variable_name + '_mann_whitney_analysis.csv')

    t_test_result = t_test.loc[t_test['pvalue'] <= 0.05]

    if len(t_test_result.index) >= 0.95 * len(t_test.index):
        print('The Mann Whitney test showed that the samples are statistically different between ' + variable_name +
              ' variable with ' + analysis_variable_name)
        t_test_result.to_csv('results/' + variable_name + '_' + analysis_variable_name +
                             '_mann_whitney_analysis_95.csv', index=False)
    else:
        print('The Mann Whitney test showed that the samples are NOT statistically different between ' +
              variable_name + ' variable with ' + analysis_variable_name)


def do_kruskal(groups, analysis_variable_values, variable_name, analysis_variable_name):
    kruskal_results = []
    for a in analysis_variable_values:
        groups_aux = [group[a] for group in groups]
        kruskal = stats.kruskal(*groups_aux)
        kruskal_results.append(kruskal)

    kruskal_test = pd.DataFrame(kruskal_results)
    kruskal_test[analysis_variable_name] = analysis_variable_values
    kruskal_test = kruskal_test[[analysis_variable_name, 'statistic', 'pvalue']].round(3)
    kruskal_test.to_csv('results/' + variable_name + '_' + analysis_variable_name + '_kruskal_analysis.csv',
                        index=False)

    # Criando gráfico
    sts = kruskal_test.columns[1:]
    qtd = len(kruskal_test.index)
    new_lines = []
    for i in range(0, qtd):
        line = kruskal_test.iloc[i]
        mw = line[0]
        stat = line[1:]

        for j in range(0, len(sts)):
            new_lines.append([mw, stat[j], sts[j]])

    kruskal_test_aux = pd.DataFrame(new_lines, columns=[analysis_variable_name, 'statistic_value', 'statistical_test'])
    kruskal_test_aux.to_csv('results/' + variable_name + '_' + analysis_variable_name + '_kruskal_analysis_aux.csv',
                            index=False)

    kruskal_data_aux = pd.read_csv('results/' + variable_name + '_' + analysis_variable_name +
                                   '_kruskal_analysis_aux.csv', header=0)

    plt.clf()
    plt.figure(figsize=(18, 8))
    kruskal_plot = sns.lineplot(x=analysis_variable_name, y="statistic_value", hue="statistical_test",
                                style="statistical_test", data=kruskal_data_aux)
    kruskal_plot.get_figure().savefig('results/' + variable_name + '_' + analysis_variable_name +
                                      '_kruskal_analysis.png')

    kruskal_test_result = kruskal_test[kruskal_test['pvalue'] <= 0.05]

    if len(kruskal_test_result.index) >= 0.95 * len(kruskal_test.index):
        print('The Kruskal-Wallis test showed that the samples are statistically different between ' + variable_name +
              ' variable with ' + analysis_variable_name)
        kruskal_test_result.to_csv('results/' + variable_name + '_' + analysis_variable_name +
                                   '_kruskal_analysis_95.csv', index=False)
    else:
        print('The Kruskal-Wallis test showed that the samples are NOT statistically different between ' +
              variable_name + ' variale with ' + analysis_variable_name)


def main():
    sns.set(style='darkgrid')
    args_parser = build_args_parser()
    args = args_parser.parse_args()

    proteins_df = pd.read_csv(args.proteins_dataset_path)
    anamnesis_df = pd.read_csv(args.anamnesis_dataset_path)
    out_path = args.out_path
    analysis_choice = args.analysis_choice
    save_graphics = args.save_graphics

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    if not os.path.exists(out_path + '/peaks'):
        os.makedirs(out_path + '/peaks')

    # Normalizando os valores das concentrações proteicas em torno de 0
    proteins_df_norm = proteins_df.loc[:, proteins_df.columns != 'sample']
    proteins_df_norm = proteins_df_norm.loc[:, proteins_df_norm.columns != 'quality']
    proteins_qtd = len(proteins_df_norm.index)
    molecular_weights = list(proteins_df_norm.columns)

    for i in range(0, proteins_qtd):
        min_value = proteins_df_norm.iloc[i].values.min()
        proteins_df_norm.iloc[i] -= min_value

    # Salvando informações dos picos
    x = proteins_df_norm.columns.values
    peaks_data = {}
    all_weights = []

    for i in range(0, proteins_qtd):
        line = proteins_df_norm.iloc[i]
        peaks, properties = find_peaks(line, prominence=10, width=5, height=25)

        if save_graphics:
            plt.figure(figsize=(12, 6))
            plt.plot(x, line)
            plt.ylabel('Fluorescence Signal')
            plt.xlabel('Protein Molecular Weight (kDa)')
            ticks, _ = plt.xticks()
            plt.xticks(ticks[::20])
            ticks, _ = plt.xticks()
            plt.xticks(ticks=ticks, labels=[j for j in x[::20]])
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
        peak_protein_weights = np.zeros(peaks.shape)

        for j in range(peaks.size):
            peak_protein_weight = float(x[peaks[j]])
            peak_value = line.values[peaks[j]]
            y_cut = line.values[cuts[j]:cuts[j + 1]]
            baseline = np.median(y_cut)
            large = np.where(y_cut > 0.5 * (peak_value + baseline))[0]
            peak_begins[j] = large.min() + cuts[j]
            peak_ends[j] = large.max() + cuts[j]
            peak_protein_weights[j] = peak_protein_weight
            areas[j] = np.sum(line.values[peak_begins[j]:peak_ends[j]] - baseline)

        all_weights.append(peak_protein_weights)

        peaks_data_partial = {}
        for k in range(0, len(peak_protein_weights)):
            peaks_data_partial[int(peak_protein_weights[k])] = [properties['peak_heights'][k],
                                                                properties['prominences'][k],
                                                                properties['width_heights'][k], areas[k]]

        peaks_data[proteins_df['sample'][i]] = peaks_data_partial

    # Preparando dados para montar planilha de picos
    all_weights = list(set([int(item) for sublist in all_weights for item in sublist]))

    peaks_database = []
    for sample in peaks_data.keys():
        curr_sample_values = []

        for weight in all_weights:
            if weight in peaks_data[sample].keys():
                sample_data = peaks_data[sample][weight]
                sample_data.insert(0, 1)
                curr_sample_values.append(sample_data)
            else:
                curr_sample_values.append([0, 0, 0, 0, 0])

        curr_sample_values = [item for sublist in curr_sample_values for item in sublist]
        curr_sample_values.insert(0, sample)
        peaks_database.append(curr_sample_values)

    peaks_cols = ["sample"]
    for weight in all_weights:
        peaks_cols.append("p_" + str(weight))
        peaks_cols.append("p_" + str(weight) + "_height")
        peaks_cols.append("p_" + str(weight) + "_prominence")
        peaks_cols.append("p_" + str(weight) + "_height")
        peaks_cols.append("p_" + str(weight) + "_area")

    # Adicionando diagnóstico a planilha de proteínas
    proteins_df_norm.insert(0, 'sample', proteins_df['sample'])
    proteins_df_norm = pd.merge(proteins_df_norm, proteins_df[['sample', 'quality']], on='sample', how='left')
    proteins_df_norm = pd.merge(proteins_df_norm, anamnesis_df[['sample', 'oral_diagnosis']], on='sample', how='left')

    # Criando a planilha com os resultados dos picos
    peaks_df = pd.DataFrame(peaks_database, columns=peaks_cols)
    peaks_df = pd.merge(peaks_df, proteins_df_norm[['sample', 'quality']], on='sample', how='left')
    peaks_df = pd.merge(peaks_df, proteins_df_norm[['sample', 'oral_diagnosis']], on='sample', how='left')

    # Salvando os resultados
    proteins_df_norm.to_csv('results/proteins_dataset_normalized.csv', index=False)
    peaks_df.to_csv('results/proteins_peaks_dataset.csv', index=False)

    proteins_df_norm_profiles = proteins_df_norm.drop(columns=['quality', 'oral_diagnosis'])
    proteins_df_norm_profiles = anamnesis_df[['sample', 'gender', 'age_group', 'oral_diagnosis']].merge(
        proteins_df_norm_profiles, on='sample', how='inner'
    )

    cols_to_drop = [c for c in peaks_df.columns if c != 'sample' and '_area' not in c]
    peaks_areas_df_profiles = peaks_df.drop(columns=cols_to_drop)
    peaks_areas_labels = peaks_areas_df_profiles.columns[1:]
    peaks_areas_labels = [label[:-5] for label in peaks_areas_labels]
    peaks_areas_labels.insert(0, 'sample')
    peaks_areas_df_profiles.columns = peaks_areas_labels
    peaks_areas_labels = peaks_areas_labels[1:]
    peaks_areas_df_profiles = anamnesis_df[['sample', 'gender', 'age_group', 'oral_diagnosis']].merge(
        peaks_areas_df_profiles, on='sample', how='inner'
    )

    # Gerando os perfis de analise por peso molecular
    profiles_lines = []
    samples_qtd = len(proteins_df_norm_profiles.index)
    for i in range(0, samples_qtd):
        line = proteins_df_norm_profiles.iloc[i]
        sample = line[0]
        gender = line[1]
        age_group = line[2]
        oral_diagnosis = line[3]
        fluorescence = line[4:]

        molecular_weights = proteins_df.columns.values[1:-1]
        molecular_weights_qtd = len(molecular_weights)
        for j in range(0, molecular_weights_qtd):
            profiles_lines.append([sample, molecular_weights[j], fluorescence[j], gender, age_group, oral_diagnosis])

    profiles_cols = ['sample', 'molecular_weight', 'fluorescence', 'gender', 'age_group', 'oral_diagnosis']
    profiles_df = pd.DataFrame(profiles_lines, columns=profiles_cols)
    profiles_df.to_csv('results/profiles_molecular_weight.csv', index=False)

    # Gerando os perfis de analise por area dos picos
    profiles_lines = []
    samples_qtd = len(peaks_areas_df_profiles.index)
    for i in range(0, samples_qtd):
        line = peaks_areas_df_profiles.iloc[i]
        sample = line[0]
        gender = line[1]
        age_group = line[2]
        oral_diagnosis = line[3]
        areas = line[4:]

        peaks = peaks_areas_df_profiles.columns.values[4:-1]
        peaks_qtd = len(peaks)
        for j in range(0, peaks_qtd):
            profiles_lines.append([sample, peaks[j], areas[j], gender, age_group, oral_diagnosis])

    profiles_cols = ['sample', 'peak', 'area', 'gender', 'age_group', 'oral_diagnosis']
    profiles_df = pd.DataFrame(profiles_lines, columns=profiles_cols)
    profiles_df.to_csv('results/profiles_peak_area.csv', index=False)

    # Cria os plots dos perfis
    profiles_df_mw = pd.read_csv('results/profiles_molecular_weight.csv', header=0)
    profiles_df_pa = pd.read_csv('results/profiles_peak_area.csv', header=0)

    if save_graphics:
        plt.clf()
        plt.figure(figsize=(18, 8))
        profiles_plot = sns.lineplot(x="molecular_weight", y="fluorescence", data=profiles_df_mw)
        profiles_plot.get_figure().savefig("results/profiles_molecular_weight.png")

        plt.clf()
        plt.figure(figsize=(18, 8))
        profiles_plot = sns.lineplot(x="peak", y="area", data=profiles_df_pa)
        profiles_plot.get_figure().savefig("results/profiles_peak_area.png")

    if analysis_choice == 'gender':
        # Análise por peso molecular
        gender_profile_df = proteins_df_norm_profiles.drop(columns=['age_group', 'oral_diagnosis', 'sample'])
        female_group = gender_profile_df[(gender_profile_df['gender'] == 'female')]
        male_group = gender_profile_df[(gender_profile_df['gender'] == 'male')]

        # Teste de Levene
        do_levene([female_group, male_group], molecular_weights, 'gender', 'molecular_weights')
        # Teste de normalidade "female"
        do_normality(female_group, molecular_weights, 'gender', 'female', 'molecular_weights')
        # Teste de normalidade "male"
        do_normality(male_group, molecular_weights, 'gender', 'male', 'molecular_weights')
        # Teste Mann Whitney
        do_mann_whitney([female_group, male_group], molecular_weights, 'gender', 'molecular_weights')

        # Análise por área de pico
        gender_profile_df = peaks_areas_df_profiles.drop(columns=['age_group', 'oral_diagnosis', 'sample'])
        female_group = gender_profile_df[(gender_profile_df['gender'] == 'female')]
        male_group = gender_profile_df[(gender_profile_df['gender'] == 'male')]

        # Teste de Levene
        do_levene([female_group, male_group], peaks_areas_labels, 'gender', 'peaks_areas')
        # Teste de normalidade "female"
        do_normality(female_group, peaks_areas_labels, 'gender', 'female', 'peaks_areas')
        # Teste de normalidade "male"
        do_normality(male_group, peaks_areas_labels, 'gender', 'male', 'peaks_areas')
        # Teste Mann Whitney
        do_mann_whitney([female_group, male_group], peaks_areas_labels, 'gender', 'peaks_areas')

        if save_graphics:
            # Plota o gráfico das florescências por gênero
            plt.clf()
            plt.figure(figsize=(18, 8))
            gender_profile_plot = sns.lineplot(x="molecular_weight", y="fluorescence", hue="gender", style="gender",
                                               data=profiles_df_mw)
            gender_profile_plot.get_figure().savefig("results/gender_profile_mw.png")

            # Plota o gráfico das áreas dos picos por gênero
            plt.clf()
            plt.figure(figsize=(18, 8))
            gender_profile_plot = sns.lineplot(x="peak", y="area", hue="gender", style="gender",
                                               data=profiles_df_pa)
            gender_profile_plot.get_figure().savefig("results/gender_profile_pa.png")
    elif analysis_choice == 'age_group':
        # Análise por peso molecular
        age_group_profile_df = proteins_df_norm_profiles.drop(columns=['gender', 'oral_diagnosis', 'sample'])
        group_13 = age_group_profile_df[(age_group_profile_df['age_group'] == '<13')]
        group_13_24 = age_group_profile_df[(age_group_profile_df['age_group'] == '13-24')]
        group_25_50 = age_group_profile_df[(age_group_profile_df['age_group'] == '25-50')]
        group_50 = age_group_profile_df[(age_group_profile_df['age_group'] == '>50')]

        # Teste de Levene
        do_levene([group_13, group_13_24, group_25_50, group_50], molecular_weights, 'age_group', 'molecular_weights')
        # Teste de normalidade "<13"
        do_normality(group_13, molecular_weights, 'age_group', '13', 'molecular_weights')
        # Teste de normalidade "13-24"
        do_normality(group_13_24, molecular_weights, 'age_group', '13_24', 'molecular_weights')
        # Teste de normalidade "25-50"
        do_normality(group_25_50, molecular_weights, 'age_group', '25_50', 'molecular_weights')
        # Teste de normalidade ">50"
        do_normality(group_50, molecular_weights, 'age_group', '50', 'molecular_weights')
        # Análise de Kruskal-Wallis
        do_kruskal([group_13, group_13_24, group_25_50, group_50], molecular_weights, 'age_group', 'molecular_weights')

        # Análise por area de pico
        age_group_profile_df = peaks_areas_df_profiles.drop(columns=['gender', 'oral_diagnosis', 'sample'])
        group_13 = age_group_profile_df[(age_group_profile_df['age_group'] == '<13')]
        group_13_24 = age_group_profile_df[(age_group_profile_df['age_group'] == '13-24')]
        group_25_50 = age_group_profile_df[(age_group_profile_df['age_group'] == '25-50')]
        group_50 = age_group_profile_df[(age_group_profile_df['age_group'] == '>50')]

        # Teste de Levene
        do_levene([group_13, group_13_24, group_25_50, group_50], peaks_areas_labels, 'age_group', 'peaks_areas')
        # Teste de normalidade "<13"
        do_normality(group_13, peaks_areas_labels, 'age_group', '13', 'peaks_areas')
        # Teste de normalidade "13-24"
        do_normality(group_13_24, peaks_areas_labels, 'age_group', '13_24', 'peaks_areas')
        # Teste de normalidade "25-50"
        do_normality(group_25_50, peaks_areas_labels, 'age_group', '25_50', 'peaks_areas')
        # Teste de normalidade ">50"
        do_normality(group_50, peaks_areas_labels, 'age_group', '50', 'peaks_areas')
        # Análise de Kruskal-Wallis
        do_kruskal([group_13, group_13_24, group_25_50, group_50], peaks_areas_labels, 'age_group', 'peaks_areas')

        if save_graphics:
            # Plota o gráfico das florescências por grupo de idade
            plt.clf()
            plt.figure(figsize=(18, 8))
            age_group_profile_plot = sns.lineplot(x="molecular_weight", y="fluorescence", hue="age_group",
                                                  style="age_group",
                                                  data=profiles_df_mw)
            age_group_profile_plot.get_figure().savefig("results/age_group_profile_mw.png")

            # Plota o gráfico das áreas dos picos por grupo de idade
            plt.clf()
            plt.figure(figsize=(18, 8))
            age_group_profile_plot = sns.lineplot(x="peak", y="area", hue="age_group",
                                                  style="age_group",
                                                  data=profiles_df_pa)
            age_group_profile_plot.get_figure().savefig("results/age_group_profile_pa.png")
    else:
        # Análise por peso molecular
        oral_diagnosis_profile_df = proteins_df_norm_profiles.drop(columns=['gender', 'age_group', 'sample'])
        group_healthy = oral_diagnosis_profile_df[(oral_diagnosis_profile_df['oral_diagnosis'] == 'healthy')]
        group_unhealthy = oral_diagnosis_profile_df[(oral_diagnosis_profile_df['oral_diagnosis'] == 'unhealthy')]

        # Teste de Levene
        do_levene([group_healthy, group_unhealthy], molecular_weights, 'oral_diagnosis', 'molecular_weights')
        # Teste de normalidade "healthy"
        do_normality(group_healthy, molecular_weights, 'oral_diagnosis', 'healthy', 'molecular_weights')
        # Teste de normalidade "unhealthy"
        do_normality(group_unhealthy, molecular_weights, 'oral_diagnosis', 'unhealthy', 'molecular_weights')
        # Análise de Kruskal-Wallis
        do_kruskal([group_healthy, group_unhealthy], molecular_weights, 'oral_diagnosis', 'molecular_weights')

        # Análise por area de pico
        oral_diagnosis_profile_df = peaks_areas_df_profiles.drop(columns=['gender', 'age_group', 'sample'])
        group_healthy = oral_diagnosis_profile_df[(oral_diagnosis_profile_df['oral_diagnosis'] == 'healthy')]
        group_unhealthy = oral_diagnosis_profile_df[(oral_diagnosis_profile_df['oral_diagnosis'] == 'unhealthy')]

        # Teste de Levene
        do_levene([group_healthy, group_unhealthy], peaks_areas_labels, 'oral_diagnosis', 'peaks_areas')
        # Teste de normalidade "healthy"
        do_normality(group_healthy, peaks_areas_labels, 'oral_diagnosis', 'healthy', 'peaks_areas')
        # Teste de normalidade "unhealthy"
        do_normality(group_unhealthy, peaks_areas_labels, 'oral_diagnosis', 'unhealthy', 'peaks_areas')
        # Análise de Kruskal-Wallis
        do_kruskal([group_healthy, group_unhealthy], peaks_areas_labels, 'oral_diagnosis', 'peaks_areas')

        if save_graphics:
            # Plota o gráfico das florescências por diaginóstico oral
            plt.clf()
            plt.figure(figsize=(18, 8))
            oral_diagnosis_profile_plot = sns.lineplot(x="molecular_weight", y="fluorescence", hue="oral_diagnosis",
                                                       style="oral_diagnosis", data=profiles_df_mw)
            oral_diagnosis_profile_plot.get_figure().savefig("results/oral_diagnosis_profile_mw.png")

            # Plota o gráfico das areas dos picos por diaginóstico oral
            plt.clf()
            plt.figure(figsize=(18, 8))
            oral_diagnosis_profile_plot = sns.lineplot(x="peak", y="area", hue="oral_diagnosis",
                                                       style="oral_diagnosis", data=profiles_df_pa)
            oral_diagnosis_profile_plot.get_figure().savefig("results/oral_diagnosis_profile_pa.png")


if __name__ == '__main__':
    main()
