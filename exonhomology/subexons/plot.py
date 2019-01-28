"""plot: Plot the MSA matrix using seaborn."""

load_plot = False
try:
    import matplotlib.pyplot as plt
    load_plot = True
except ImportError as err:
    print("matplotlib.pyplot loading has failed with ImportError: ", err)
    if 'tkinter' in str(err):
        print("""
            You should install tkinter for Python 3 to access the
            subexons.plot.plot_msa_subexons function.
            In Ubuntu/Debian you can do: sudo apt-get install python3-tk
            """)

if load_plot:
    import numpy as np
    import pandas as pd
    import seaborn

    # import matplotlib.pyplot as plt gives the following ImportError:
    # No module named '_tkinter', please install the python3-tk package


    def plot_msa_subexons(gene_ids,
                          msa_matrix,
                          subexon_table=None,
                          outfile='alignment.png'):
        """
        Save a plot from the msa matrix where cells are subexons.

        If `subexon_table` is `None`, the binary plot with the constitutive
        exons is not generated. Otherwise, it'll be in the `_constitutive.png`
        file. The `subexon_table` also allows the coloring by Subexon ID for
        the `_subexon_id.png` file.
        """
        # stack... stop-seaborn-plotting-multiple-figures-on-top-of-one-another
        plt.figure()

        subexon_number = int(msa_matrix[~np.isnan(msa_matrix)].max()) + 1
        cmap = seaborn.color_palette('tab10', subexon_number)

        msa_matrix_df = pd.DataFrame(msa_matrix)
        msa_matrix_df['Gene ID'] = gene_ids
        msa_matrix_df.set_index('Gene ID', inplace=True)
        sns_plot = seaborn.heatmap(msa_matrix_df, cmap=cmap)

        # This depends on the seaborn version: seaborn-0.8.1
        sns_plot.get_figure().savefig(outfile, bbox_inches='tight')
        # stackoverflow: second-y-axis-label-getting-cut-off

        if subexon_table is not None:
            msa_matrix_df_copy = msa_matrix_df.copy()

            id2fraction = pd.Series(
                subexon_table['Transcript fraction'].values,
                index=subexon_table['Subexon Index']).to_dict()

            subexon_subdf = subexon_table.loc[:, [
                'Number of transcripts for subexon', 'Transcripts in gene',
                'Subexon Index', 'Subexon ID'
            ]]

            subexon_subdf = subexon_subdf.drop_duplicates().set_index(
                'Subexon Index')

            nrow, ncol = msa_matrix_df.shape

            for i in range(0, nrow):
                for j in range(0, ncol):
                    subexon_index = msa_matrix_df_copy.iloc[i, j]
                    if not np.isnan(subexon_index):
                        subexon_index = int(subexon_index)
                        fraction = id2fraction[subexon_index]
                        if fraction == 1.0:
                            msa_matrix_df_copy.iloc[i, j] = 1.0
                        else:
                            nt = subexon_subdf.loc[subexon_index][
                                'Number of transcripts for subexon']
                            n = subexon_subdf.loc[subexon_index][
                                'Transcripts in gene']
                            if n == 1:
                                msa_matrix_df_copy.iloc[i, j] = 0.33
                            elif n - 1 == nt:
                                msa_matrix_df_copy.iloc[i, j] = 0.66
                            else:
                                msa_matrix_df_copy.iloc[i, j] = 0.0

            plt.figure()
            outfile_1 = outfile.replace('.png', '_constitutive.png')
            sns_const_plot1 = seaborn.heatmap(
                msa_matrix_df_copy, cmap="RdYlGn")
            sns_const_plot1.get_figure().savefig(
                outfile_1, bbox_inches='tight')

            msa_matrix_df_copy = msa_matrix_df.copy()
            for i in range(0, nrow):
                for j in range(0, ncol):
                    subexon_index = msa_matrix_df_copy.iloc[i, j]
                    if not np.isnan(subexon_index):
                        subexon_index = int(subexon_index)
                        msa_matrix_df_copy.iloc[i, j] = subexon_subdf.loc[
                            subexon_index]['Subexon ID']
            plt.figure()
            outfile2 = outfile.replace('.png', '_subexon_id.png')
            sns_const_plot2 = seaborn.heatmap(
                msa_matrix_df_copy, cmap="Spectral")
            sns_const_plot2.get_figure().savefig(outfile2, bbox_inches='tight')

        return sns_plot
