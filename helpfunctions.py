import sys
import openpyxl
import itertools
import numpy as np
import pandas as pd
import validation as valid
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication
from groupSeperator import AssignValuesWindow


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
    print(f"The SVG file '{title}' saved successfully")
    plt.close()


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

    control_types = ['CONTROL', 'DMSO', 'PBS']
    compound_name = pairs_df['compound_name'].unique()

    control_list = list(set(compound_name) & set(control_types))
    inhibitor_list = list(set(compound_name) - set(control_list))

    app = QApplication(sys.argv)
    window = AssignValuesWindow(control_list, inhibitor_list, cell_name)
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
                worksheet.freeze_panes = "A2"
            print(f"The sheet '{sheet_name}' created successfully")
        except Exception as e:
            print(f"Error occurred while creating the sheet: {e}")
    else:
        print(f"The sheet '{sheet_name}' was not created because the DataFrame is empty")
