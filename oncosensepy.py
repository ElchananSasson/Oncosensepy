import os
import sys
import numpy as np
import pandas as pd
import exceptions as e
import helpfunctions as hf
import validation as valid
from scipy.stats import ttest_ind
from openpyxl import load_workbook
from PyQt5.QtWidgets import QApplication
from cellNamesGUI import AssignNamesValuesWindow


def get_LGE_data(data_set_path):
    """
    The function reads Excel sheets ('L' and 'G') from the specified file path and returns clear DataFrames without
    missing values.

    Param:
        data_set_path (str): The path to the Excel file containing the data.

    Returns:
        l_df (pandas.DataFrame): A DataFrame containing the data from the 'L' sheet.
        g_df (pandas.DataFrame): A DataFrame containing the data from the 'G' sheet.
    """
    l_df = pd.read_excel(data_set_path, sheet_name='L').fillna(0)
    g_df = pd.read_excel(data_set_path, sheet_name='G').fillna(0)
    err_limit_lambda = pd.read_excel(data_set_path, sheet_name='ErrorLimitLambda').columns.values[0]

    l_df['compound_name'] = l_df['compound_name'].apply(lambda x: 'CONTROL' if x == 0 else x)
    l_df['2D_3D'] = l_df['2D_3D'].apply(lambda x: '-0-' if x == 0 else x)
    l_df['dosage'] = l_df['dosage'].apply(lambda x: '-0-' if x == 0 else x)
    l_df['time'] = l_df['time'].apply(lambda x: '0hr' if x == 0 else x)
    if 0 in l_df['cell_line_name'].values:
        e.InvalidCellLineException("Cell line name has missing values")

    if 0 in g_df['UID'].values:
        e.InvalidUIDException("UID has missing values")

    return l_df, g_df, err_limit_lambda


def important_L(l_df, err_limit, threshold, new_sheet=False, sheet_name='important_L', data_path=''):
    """
    This method returns a DataFrame with only the important columns. An important column is determined
    by whether the number of cells whose value is higher in absolute value than the error limit,
    is greater than or equal to the threshold.

    Params:
        l_df (pandas.DataFrame): The DataFrame to be checked.
        err_limit (float): The error limit.
        threshold (int): The number of significant values.
        new_sheet (bool): If True, creates a new sheet. Default is False.
        sheet_name (str): The name of the sheet to be created. Default is 'important_L'.
        path (str): The path where the new sheet will be created. Default is an empty string.

    Returns:
        pandas.DataFrame: The DataFrame with only the important columns selected.
    """
    valid.is_valid_L(l_df)
    if threshold < 0:
        raise e.NegativeNumberException("Threshold should be positive number")
    new_df = l_df.loc[:, :5].copy()
    for i in range(1, len(l_df.columns) - 5):
        count = 0
        for cell in l_df[i]:
            if abs(cell) > err_limit:
                count += 1
        if count >= threshold:
            new_df[i] = l_df[i].copy()

    if new_sheet:
        print(f"Creating '{sheet_name}'..")
        hf.create_new_sheet(new_df, data_path, sheet_name)

    return new_df


def filter_by_col(df, col, filter_list, new_sheet=False, sheet_name='filter_by_col', data_path=''):
    """
    This function filters data by a certain column and by a list of values it receives.

    Params:
        df (pandas.DataFrame): The DataFrame to filter.
        col (str): The column according to which the filtering will be performed.
        filter_list (list): A list of values that we would like to appear in the selected column.
        new_sheet (bool): If True, creates a new sheet. Default is False.
        sheet_name (str): The name of the sheet to be created. Default is 'filter_by_col'.
        path (str): The path where the new sheet will be created. Default is an empty string.

    Returns:
        pandas.DataFrame: The DataFrame after filtering.
    """
    valid.is_valid_L(df)
    filter_df = df.loc[(df[col].isin(filter_list))]
    if len(filter_df) == 0:
        print(f"There is no data to show by '{col}' filtering")

    if new_sheet:
        print(f"Creating '{sheet_name}'..")
        hf.create_new_sheet(filter_df, data_path, sheet_name)

    return filter_df


def sort_G_values(g_df, cols, save_plot_path='', save=False, new_sheet=False, sheet_name='Sorted_G', data_path=''):
    """
    This function accepts columns representing processes and sorts for each process its proteins.
    In addition, the function saves the plot of each process if the 'save' argument is set to True.
    Otherwise, the plots will be displayed one by one.

    Params:
        g_df (pandas.DataFrame): The G_values DataFrame.
        cols (list): The process to sort its G_values.
        save_plot_path (str): The path to save the figures. If not specified or set to an empty string (''),
                              the function will not save the plots.
        save (bool): If True, the function will save the plot of each process in the specified path.
        new_sheet (bool): If True, creates a new sheet. Default is False.
        sheet_name (str): The name of the sheet to be created. Default is 'Sorted_G'.
        data_path (str): The path where the new sheet will be created. Default is an empty string.

    Returns:
        pandas.DataFrame: Sorted G_values.
    """

    if save_plot_path == '':
        save_plot_path = os.getcwd()
    valid.is_valid_path(save_plot_path)
    important_g = pd.DataFrame(columns=pd.MultiIndex.from_product([cols, ['UID', 'Effect']]))
    names_g, values_g = {}, {}
    for col in important_g.columns:
        if col[1] == 'UID':
            names_g = g_df[col[1]].to_dict()
        else:
            values_g = g_df[col[0]].to_dict()
            sorted_values_g = dict(sorted(values_g.items(), key=lambda item: item[1]))
            list_values_g = list(sorted_values_g.values())

            sorted_names_g = dict(sorted(names_g.items(), key=lambda item: sorted_values_g[item[0]]))
            list_names_g = list(sorted_names_g.values())

            if save:
                hf.plot_G_values(f'Process {col[0]}', list_names_g, list_values_g, save_plot_path)

            important_g[(col[0], 'UID')] = list_names_g
            important_g[col] = list_values_g

    if new_sheet:
        print(f"Creating '{sheet_name}'..")
        hf.create_new_sheet(important_g, data_path, sheet_name)

    return important_g


def analyze_pairs(important_l, cell_line_list=None, fixed_col='time', p_value=0.05, only_avg=False, data_path=''):
    """
    This function analyzes pairs of compounds in a dictionary of Pandas dataframes.

    Params:
        important_l (pandas.DataFrame): The DataFrame with only the important columns.
        cell_line_list (list): A list of cell lines to analyze.
        fixed_col (str): The name of the column that will remain fixed in each pair. Default is 'time'.
        p_value (float): The p-value threshold for determining whether the difference
                         between means is significant. Default is 0.05.
        only_avg (bool): If True, return only the averages for each column. Default is False.
        data_path (str): The path where the updated dataframes will be saved. Default is an empty string.
    """

    workbook = load_workbook(filename=data_path)
    sheet = workbook['ErrorLimitLambda']
    error_value = sheet['A1'].value
    if cell_line_list is None:
        cell_line_list = important_l['cell_line_name'].unique().tolist()
        app_names = QApplication(sys.argv)
        name_window = AssignNamesValuesWindow(cell_line_list)
        name_window.show()
        app_names.exec_()
        cell_line_list = name_window.result
        del app_names

    for cell_line in cell_line_list:
        print(f"analyzing '{cell_line}' cell line..")
        pairs_dict = hf.pairs_df_to_dict(important_l, cell_line, fixed_col=fixed_col)

        keys_to_remove, compound_names = [], []
        averages = {}

        for key, sub_df in pairs_dict.items():
            col_names = sub_df.columns.tolist()
            time_col_idx = col_names.index('time')
            analysis_cols = [c for c in col_names[time_col_idx + 1:]]

            dfs_to_concat = []

            for col in analysis_cols:
                if key[1] == key[2]:
                    df_first = sub_df.loc[sub_df[fixed_col] == key[3], col]
                    df_second = sub_df.loc[sub_df[fixed_col] == key[4], col]
                else:
                    df_first = sub_df.loc[sub_df['compound_name'] == key[1], col]
                    df_second = sub_df.loc[sub_df['compound_name'] == key[2], col]

                sign_changed = False
                df_first_mean = df_first.mean()
                df_second_mean = df_second.mean()
                if np.sign(df_first_mean) != np.sign(df_second_mean):
                    sign_changed = True
                if df_first.shape[0] == 1 or df_second.shape[0] == 1:
                    p = p_value + 1
                    # TODO schedule meeting with Nissim about it (the problem is that we have 1 value)
                else:
                    t, p = ttest_ind(df_first, df_second)
                if sign_changed or (p <= p_value):
                    if (abs(df_first_mean) > error_value) or (abs(df_second_mean) > error_value):
                        averages[col] = (df_first_mean, df_second_mean)
                        add_res = pd.DataFrame([[hf.add_reason(sign_changed, p, p_value)]], columns=[col])
                        add_res = add_res.rename(index={0: 'Reason'})
                        df_with_res_row = pd.concat([sub_df[[col]], add_res])
                        dfs_to_concat.append(df_with_res_row)

            if len(dfs_to_concat) == 0:
                keys_to_remove.append(key)
            else:
                if only_avg:
                    updated_dfs = []
                    for df_col in dfs_to_concat:
                        last_row = df_col.iloc[-1].values[0]
                        compound_names = sub_df['compound_name'].dropna().unique().tolist()
                        if len(compound_names) == 1:
                            compound_names.append(compound_names[0])
                        new_col = pd.DataFrame(
                            [averages[df_col.columns[0]][0], averages[df_col.columns[0]][1], last_row],
                            columns=[df_col.columns[0]],
                            index=[compound_names[0] + " AVG", compound_names[1] + " AVG", "Reason"])

                        updated_dfs.append(new_col)
                    dfs_to_concat = updated_dfs

                new_df = pd.concat(dfs_to_concat, axis=1)
                new_df = new_df.reindex(sorted(new_df.columns), axis=1)

                if only_avg:
                    columns_to_select = ['cell_line_name', 'compound_name', fixed_col]

                    index_first_time = sub_df[fixed_col].str.extract(r'(\d+)', expand=False).astype(int).idxmin()
                    index_second_time = sub_df[fixed_col].str.extract(r'(\d+)', expand=False).astype(int).idxmax()

                    df_without_barcode = pd.DataFrame(columns=columns_to_select)
                    df_without_barcode = pd.concat(
                        [df_without_barcode] + [pd.DataFrame(columns=df_without_barcode.columns)] * 3,
                        ignore_index=True)

                    df_without_barcode.loc[0] = sub_df.loc[index_first_time, columns_to_select]
                    df_without_barcode.loc[1] = sub_df.loc[index_second_time, columns_to_select]

                    if df_without_barcode.shape[0] == 2:
                        last_row = df_without_barcode.iloc[-1]
                        new_row = pd.DataFrame([last_row], columns=df_without_barcode.columns)
                        df_without_barcode = pd.concat([df_without_barcode, new_row], ignore_index=True)
                    compound_names.append('')

                    df_without_barcode['compound_name'] = compound_names
                    df_without_barcode.index = [compound_names[0] + " AVG", compound_names[1] + " AVG", "Reason"]
                    df_without_barcode.loc["Reason"] = np.nan
                    pairs_dict[key] = pd.concat([df_without_barcode, new_df], axis=1)
                else:
                    pairs_dict[key] = pd.concat(
                        [sub_df[['cell_line_name', 'compound_name', '2D_3D', 'dosage', 'time']], new_df],
                        axis=1)

        for key in keys_to_remove:
            pairs_dict.pop(key)

        pairs_df = hf.create_pairs_df(pairs_dict)
        sheet_name = cell_line + "_FULL"
        if only_avg:
            sheet_name = cell_line + "_avg_by_" + fixed_col

        hf.create_new_sheet(pairs_df, data_path, sheet_name)


def analyze_control_treatment(df, cell_name, control_list=None, p_value=0.05):
    """
    This function analyzes a single Pandas dataframe for a specified cell line and performs a t-test between the
    means of the control and treatment groups for each time-point column. If the difference between means is not
    significant (as determined by the p-value threshold), the column is dropped from the dataframe. The default
    control compounds are 'CONTROL', 'DMSO', and 'PBS'.

    Params:
    df (Pandas dataframe): The input dataframe to be analyzed.
    cell_name (str): The name of the cell line to be analyzed.
    control_list (list of str): A list of control compound names. Default is ['CONTROL', 'DMSO', 'PBS'].
    p_value (float): The p-value threshold for determining whether the difference between means is significant.
    Default is 0.05.

    Returns:
    Pandas dataframe: A copy of the input dataframe with any columns where the difference between means was not
    significant (as determined by the p-value threshold) dropped.
    """
    if control_list is None:
        control_list = ['CONTROL', 'DMSO', 'PBS']

    new_df = df.loc[df['cell_line_name'] == cell_name]

    col_names = new_df.columns.tolist()
    time_col_idx = col_names.index('time')
    analysis_cols = [c for c in col_names[time_col_idx + 1:]]

    for col in analysis_cols:
        df_first = new_df.loc[new_df['compound_name'].isin(control_list), col]
        df_second = new_df.loc[~new_df['compound_name'].isin(control_list), col]

        first_avg = df_first.mean()
        second_avg = df_second.mean()

        t, p = ttest_ind(df_first, df_second)
        if (np.sign(first_avg) == np.sign(second_avg)) and (p > p_value):
            new_df = new_df.drop(col, axis=1)

    return new_df
