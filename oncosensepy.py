import os
import sys
import openpyxl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from PyQt5.QtWidgets import QApplication
import validation as valid
from scipy.stats import ttest_ind
import exceptions as e
from groupSeperator import AssignValuesWindow


def get_LG_data(data_set_path):
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

    l_df['compound_name'] = l_df['compound_name'].apply(lambda x: 'CONTROL' if x == 0 else x)
    l_df['time'] = l_df['time'].apply(lambda x: '0hr' if x == 0 else x)
    if 0 in l_df['cell_line_name'].values:
        e.InvalidCellLineException("Cell line name has missing values")

    if 0 in g_df['UID'].values:
        e.InvalidUIDException("UID has missing values")

    return l_df, g_df


def important_L(df, err_limit, threshold):
    """
        This method returns a DataFrame with only the important columns. An important column is determined
        by whether the number of cells whose value is higher in absolute value than the error limit,
        is greater than or equal to the threshold.

        Params:
            df (pandas.DataFrame): The DataFrame to be checked.
            err_limit (float): The error limit.
            threshold (int): The number of significant values.

        Returns:
            pandas.DataFrame: The DataFrame with only the important columns selected.
    """
    valid.is_valid_L(df)
    if threshold < 0:
        raise e.NegativeNumberException("Threshold should be positive number")
    new_df = df.loc[:, :5].copy()
    for i in range(1, len(df.columns) - 5):
        count = 0
        for cell in df[i]:
            if abs(cell) > err_limit:
                count += 1
        if count >= threshold:
            new_df[i] = df[i].copy()
    return new_df


def filter_by_col(df, col, filter_list):
    """
        This function filters data by a certain column and by a list of values it receives.

        Params:
            df (pandas.DataFrame): The DataFrame to filter.
            col (str): The column according to which the filtering will be performed.
            filter_list (List): A list of values that we would like to appear in the selected column.

        Returns:
            pandas.DataFrame: The DataFrame after filtering.
    """
    valid.is_valid_L(df)
    filter_df = df.loc[(df[col].isin(filter_list))]
    if len(filter_df) == 0:
        print(f"There is no data to show by '{col}' filtering")
    return filter_df


def find_edges(list_names_g, list_values_g):
    """
       The function finds the start and end edges of a graph represented as a DataFrame.

       Params:
           list_values_g(List[float]): The list of the values
           list_names_g(List[str]): The list of the names of the proteins

       Returns:
           Tuple[List[str], List[str]]: A tuple containing the start and end edges of the graph as lists of strings.
    """
    if len(list_values_g) == 0:
        return [], []
    # split the list into 2 lists
    half_ind = len(list_values_g) // 2
    first_half = list_values_g[:half_ind]
    second_half = list_values_g[half_ind:]
    distance, max_distance = 0, 0
    first_indices, second_indices = [], []

    # find the max destination between 2 points and save their index and the value in first_indices
    for i in range(1, len(first_half)):
        distance = abs(abs(first_half[i]) - abs(first_half[i - 1]))
        if distance > max_distance:
            max_distance = distance
            first_indices = [i - 1, i, max_distance]

    # find the max destination between 2 points and save their index and the value in second_indices
    max_distance = 0
    for i in range(1, len(second_half)):
        distance = abs(abs(second_half[i]) - abs(second_half[i - 1]))
        if distance > max_distance:
            max_distance = distance
            second_indices = [i - 1, i, max_distance]

    lower_edge = list_names_g[: first_indices[1]]
    upper_edge = list_names_g[half_ind + second_indices[1]:]
    return lower_edge, upper_edge


def sort_G_values(g_df, cols, path='', save=False):
    """
        This function accepts columns representing processes and sorts for each process its proteins.
        In addition, the function saves the plot of each process if the 'save' argument is set to True.
        Otherwise, the plots will be displayed one by one.

        Params:
            g_df (pandas.DataFrame): The G_values DataFrame.
            cols (list): The process to sort its G_values.
            path (str): The path to save the figures. If not specified or set to an empty string (''),
                        the function will not save the plots.
            save (bool): If True, the function will save the plot of each process in the specified path.

        Returns:
            pandas.DataFrame: Sorted G_values.
    """

    def plot_G_values(title, uid, values, save_path):
        """
            This function accepts columns representing processes and sorts for each process its proteins.
            In addition, the function saves the plot of process

            Params:
                title (str): The plot title.
                uid (list): The sorted list of G_UID.
                values (list): The sorted list of G_values.
                save_path (str): The path to save the figures, if None the plots will be displayed one by one
        """
        valid.is_valid_path(save_path)

        plt.title(title)
        plt.figure(figsize=(50, 30))
        plt.scatter(uid, values)
        plt.xticks(uid, [f'{name} ({i})' for i, name in enumerate(uid)], rotation=90, fontsize=6)
        plt.ylabel('Effect')

        for i, txt in enumerate(range(len(uid))):
            plt.text(uid[i], values[i] + 0.002, str(i), fontsize=5)

        lower_edge, upper_edge = find_edges(uid, values)
        if lower_edge != [] and upper_edge != []:
            lower_edge_indices = [uid.index(val) for val in lower_edge]
            upper_edge_indices = [uid.index(val) for val in upper_edge]
            plt.scatter([uid[i] for i in lower_edge_indices], [values[i] for i in lower_edge_indices], color='red')
            plt.scatter([uid[i] for i in upper_edge_indices], [values[i] for i in upper_edge_indices], color='red')

        plt.savefig(save_path + '/' + f'{title}.SVG', dpi=300)
        plt.close()

    # plot_G_values function
    if path == '':
        path = os.getcwd()
    valid.is_valid_path(path)
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
                plot_G_values(f'Process {col[0]}', list_names_g, list_values_g, path)

            important_g[(col[0], 'UID')] = list_names_g
            important_g[col] = list_values_g

    return important_g


def pairs_df_to_dict(df, cell_name, control_list=None, inhibitor_list=None, fixed_col='time'):
    """
    Convert a Pandas dataframe to a dictionary of pairs of dataframes.

    Params:
        df (pandas.DataFrame): The input dataframe to convert.
        cell_name (str): The name of the cell line to filter the dataframe by.
        control_list (list): A list of control compound names.
        inhibitor_list (list): A list of inhibitor compound names.
        fixed_col (str): The name of the column that will remain fixed in each pair. Default is 'time'.

    Returns:
        dict: A dictionary where each key is a tuple representing a pair of compounds, along with optional time points.
              The corresponding value is a Pandas dataframe containing data for that pair.
              The keys are generated by combining cell_name, compound names and time points.
              There are two types of keys:
              (1) Keys for pairs of control_list and inhibitor_list with same fixed_col.
              (2) Keys for pairs of inhibitor_list with itself with different fixed_col.
    """
    # Filter the input dataframe to only include rows with the specified cell name
    pairs_df = df.loc[df['cell_line_name'] == cell_name]

    if control_list is None or inhibitor_list is None:
        if control_list is None:
            control_list = ['CONTROL', 'DMSO', 'PBS']

        if inhibitor_list is None:
            compound_name = pairs_df['compound_name'].unique()
            inhibitor_list = list(set(compound_name) - set(control_list))
        app = QApplication(sys.argv)
        window = AssignValuesWindow(control_list, inhibitor_list)
        window.show()
        app.exec_()
        control_list, inhibitor_list = window.result

    full_list = control_list + inhibitor_list

    # Filter the dataframe to only include rows with compound names in the full list
    pairs_df = pairs_df.loc[pairs_df['compound_name'].isin(full_list)]

    pairs_dict = {}
    # Pairs of control_list and inhibitor_list with same fixed_col
    for i, j in itertools.product(control_list, inhibitor_list):
        for col in pairs_df[fixed_col].unique():
            df_i_j_t = pairs_df.loc[(pairs_df['compound_name'] == i) & (pairs_df[fixed_col] == col)]
            df_j_i_t = pairs_df.loc[(pairs_df['compound_name'] == j) & (pairs_df[fixed_col] == col)]
            if not (df_i_j_t.empty or df_j_i_t.empty):
                pairs_dict[(cell_name, i, j, col)] = pd.concat([df_i_j_t, df_j_i_t])

    # Pairs of inhibitor_list with itself with different fixed_col
    for i in inhibitor_list:
        unique_fixed_col_i = pairs_df.loc[pairs_df['compound_name'] == i, fixed_col].unique()
        for t1, t2 in itertools.combinations(unique_fixed_col_i, 2):
            df_i_t1 = pairs_df.loc[(pairs_df['compound_name'] == i) & (pairs_df[fixed_col] == t1)]
            df_i_t2 = pairs_df.loc[(pairs_df['compound_name'] == i) & (pairs_df[fixed_col] == t2)]
            if not (df_i_t1.empty and df_i_t2.empty):
                pairs_dict[(cell_name, i, i, t1, t2)] = pd.concat([df_i_t1, df_i_t2])

    return pairs_dict


def analyze_pairs(pairs_dict, p_value=0.05, fixed_col='time', display=False, only_avg=False):
    """
    This function analyzes pairs of compounds in a dictionary of Pandas dataframes.

    Params:
        pairs_dict (dict): A dictionary containing Pandas dataframes for each pair of compounds.
                           The keys of the dictionary are tuples of two compound names.
        p_value (float): The p-value threshold for determining whether the difference
                         between means is significant. Default is 0.05.
        fixed_col (str): The name of the column that will remain fixed in each pair. Default is 'time'.
        display (bool): If True, the function print the dictionary
        only_avg (bool): If True, return only the averages for each col

    Returns:
        dict: A dictionary with the same keys as the input dictionary, but with updated dataframes.
              If a key-value pair is removed from the dictionary because none of the columns pass the test,
              it will not appear in the output dictionary.
    """
    keys_to_remove, compound_names = [], []
    averages = {}
    for key, df in pairs_dict.items():
        col_names = df.columns.tolist()
        time_col_idx = col_names.index('time')
        analysis_cols = [c for c in col_names[time_col_idx + 1:]]

        dfs_to_concat = []
        for col in analysis_cols:
            if key[1] == key[2]:
                df_first = df.loc[df[fixed_col] == key[3], col]
                df_second = df.loc[df[fixed_col] == key[4], col]
            else:
                df_first = df.loc[df['compound_name'] == key[1], col]
                df_second = df.loc[df['compound_name'] == key[2], col]

            sign_changed = False
            if np.sign(df_first.mean()) != np.sign(df_second.mean()):
                sign_changed = True

            t, p = ttest_ind(df_first, df_second)
            if sign_changed or (p <= p_value):
                averages[col] = (df_first.mean(), df_second.mean())

                # add reason row
                add_res = pd.DataFrame([[add_reason(sign_changed, p, p_value)]], columns=[col])
                add_res = add_res.rename(index={0: 'Reason'})

                # append the new row to df[[col]]
                df_with_res_row = pd.concat([df[[col]], add_res])
                dfs_to_concat.append(df_with_res_row)

        if len(dfs_to_concat) == 0:
            keys_to_remove.append(key)
        else:
            if only_avg:
                updated_dfs = []
                for df_col in dfs_to_concat:
                    last_row = df_col.iloc[-1].values[0]
                    compound_names = df['compound_name'].dropna().unique().tolist()
                    if len(compound_names) == 1:
                        compound_names.append(compound_names[0])
                    new_col = pd.DataFrame([averages[df_col.columns[0]][0], averages[df_col.columns[0]][1], last_row],
                                           columns=[df_col.columns[0]], index=[compound_names[0] + " AVG",
                                                                               compound_names[1] + " AVG", "Reason"])

                    updated_dfs.append(new_col)
                dfs_to_concat = updated_dfs
            new_df = pd.concat(dfs_to_concat, axis=1)
            new_df = new_df.reindex(sorted(new_df.columns), axis=1)
            if only_avg:
                columns_to_select = ['cell_line_name', 'compound_name', '2D_3D', 'dosage', 'time']
                df_without_barcode = df.loc[df.index[:3], columns_to_select]
                compound_names.append('')

                df_without_barcode['compound_name'] = compound_names
                df_without_barcode.index = [compound_names[0] + " AVG", compound_names[1] + " AVG", "Reason"]
                df_without_barcode.loc["Reason"] = np.nan
                pairs_dict[key] = pd.concat([df_without_barcode, new_df], axis=1)
            else:
                pairs_dict[key] = pd.concat(
                    [df[['barcode', 'cell_line_name', 'compound_name', '2D_3D', 'dosage', 'time']], new_df], axis=1)

    for key in keys_to_remove:
        pairs_dict.pop(key)

    if display:
        pd.set_option("display.max_rows", None)  # Display all rows
        pd.set_option("display.max_columns", None)  # Display all columns
        print(pairs_dict)

    return pairs_dict


def add_reason(sign_changed, p, p_value):
    """
        Adds a reason row to a Pandas DataFrame indicating the result of the analysis.

        Params:
            df (DataFrame): The input DataFrame.
            col (str): The name of the column in the DataFrame to add the reason row.
            total_avg (float): The average value calculated for the column across groups.
            p (float): The p-value calculated for the statistical test.
            p_value (float): The p-value threshold for determining significance.
            col_len (int): The length of the column.

        Returns:
            pandas.DataFrame: The input DataFrame with the reason row added.
        """
    if sign_changed and (p <= p_value):
        return "P-Value and Sign change"

    elif sign_changed:
        return "Sign change"

    elif p <= p_value:
        return "P-Value"


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

    # Filter the input dataframe to only include rows with the specified cell name
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


def create_pairs_df(pairs_dict):
    """
    This function takes a dictionary of paired DataFrames and concatenates them into a single DataFrame for comparison.

    Param:
        pairs_dict (dict): A dictionary of paired DataFrames.

    Returns:
        pandas.DataFrame: A concatenated DataFrame of paired DataFrames for comparison.
    """
    comp_list = []
    for i, df in enumerate(pairs_dict.values()):
        comp_list.append(df)
        if i < len(pairs_dict) - 1:
            comp_list.append(pd.DataFrame(np.nan, index=['-'], columns=df.columns))

    pairs_df = pd.concat(comp_list, sort=False)
    return pairs_df


def create_csv(df, name='default_name', path=None):
    """
        This method inserts the DataFrame into a csv file.

        Params:
            df (pandas.DataFrame): The DataFrame to insert into a new csv.
            name (str): File name.
            path (str): The file path.
    """
    if df.empty:
        print(f"The csv file '{name}' was not created because the DataFrame is empty")
        return

    # Create the folder if it doesn't exist
    elif path is None:
        folder = 'csv'
        if not os.path.exists(folder):
            os.makedirs(folder)
        path = folder + '/' + name + '.csv'  # Save the data frame to the CSV file
        df.to_csv(path)
    else:
        df.to_csv(path + '/' + name + '.csv')
    print(f"The csv '{name}' created successfully")


def create_new_sheet(df, path, sheet_name):
    """
        This method inserts the DataFrame into a new sheet in an existing Excel file.

        Params:
            df (pandas.DataFrame): The DataFrame to insert into a new sheet.
            path (str): The file path.
            sheet_name (str): The name of the new sheet.
    """
    if not df.empty:
        try:
            wb = openpyxl.load_workbook(path)
            if sheet_name in wb.sheetnames:
                sheet = wb[sheet_name]
                wb.remove(sheet)
                wb.save(path)
        except Exception as e:
            print(f"Error occurred while removing the sheet: {e}")
            return
        try:
            with pd.ExcelWriter(path, mode='a') as writer:
                df.to_excel(writer, sheet_name=sheet_name)
                workbook = writer.book
                worksheet = workbook[sheet_name]

                # Freeze the top row
                worksheet.freeze_panes = "A2"
            print(f"The sheet '{sheet_name}' created successfully")
        except Exception as e:
            print(f"Error occurred while creating the sheet: {e}")
    else:
        print(f"The sheet '{sheet_name}' was not created because the DataFrame is empty")
