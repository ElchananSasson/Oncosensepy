import itertools
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import validation as valid


def important_L(df, err_limit, threshold):
    """
        This method returns a DataFrame with only the important columns. An important column is determined
        by whether the number of cells whose value is higher in absolute value than the error limit,
        is greater than or equal to the threshold.

        Arguments:
            df (pandas.DataFrame): The DataFrame to be checked.
            err_limit (float): The error limit.
            threshold (int): The number of significant values.

        Returns:
            pandas.DataFrame: The DataFrame with only the important columns selected.
    """
    valid.is_valid_L(df)
    new_df = df.loc[:, :5].copy()
    new_df['time'] = new_df['time'].apply(lambda x: '0hr' if x == 0 else x)
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

        Arguments:
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


def plot_G_values(title, uid, values, path):
    """
        This function accepts columns representing processes and sorts for each process its proteins.
        In addition, the function saves the plot of process

        Arguments:
            title (str): The plot title.
            uid (list): The list of G_UID.
            values (list): The list of G_values.
            path (str): The path to save the figures, if None the plots will be displayed one by one
    """
    valid.is_valid_path(path)
    plt.title(title)
    plt.figure(figsize=(50, 30))
    plt.scatter(uid, values)
    plt.xticks(uid, [f'{name} ({i})' for i, name in enumerate(uid)], rotation=90, fontsize=6)
    plt.ylabel('Effect')

    for i, txt in enumerate(range(len(uid))):
        plt.text(uid[i], values[i] + 0.002, str(i), fontsize=5)

    plt.savefig(path + '/' + f'{title}.SVG', dpi=300)
    plt.close()


def sort_plot_G_values(g_df, cols, path=''):
    """
        This function accepts columns representing processes and sorts for each process its proteins.
        In addition, the function saves the plot of each process

        Arguments:
            g_df (pandas.DataFrame): The G_values DataFrame.
            cols (list): The process to sort it's G_values.
            path (str): The path to save the figures, if None the plots will be displayed one by one

        Returns:
            pandas.DataFrame: Sorted G_values.
    """
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

            plot_G_values(f'Process {col[0]}', list_names_g, list_values_g, path)

            important_g[(col[0], 'UID')] = list_names_g
            important_g[col] = list_values_g

    return important_g


def create_pairs_df_dict(df, cell_name, control_list, inhibitor_list, fixed_col='time'):
    # Filter the input dataframe to only include rows with the specified cell name
    df = df.loc[df['cell_line_name'] == cell_name]

    full_list = control_list + inhibitor_list

    # Filter the dataframe to only include rows with compound names in the full list
    df = df.loc[df['compound_name'].isin(full_list)]

    # If the filtered dataframe is empty, print an error message and exit the program
    if df.empty:
        print(f'There is a problem here: {control_list} or here: {inhibitor_list}')
        sys.exit()

    pairs_dict = {}
    # Pairs of control_list and inhibitor_list with same fixed_col
    for i, j in itertools.product(control_list, inhibitor_list):
        for col in df[fixed_col].unique():
            df_i_j_t = df.loc[(df['compound_name'] == i) & (df[fixed_col] == col)]
            df_j_i_t = df.loc[(df['compound_name'] == j) & (df[fixed_col] == col)]
            if not (df_i_j_t.empty or df_j_i_t.empty):
                pairs_dict[(cell_name, i, j, col)] = pd.concat([df_i_j_t, df_j_i_t])

    # Pairs of inhibitor_list with itself with different fixed_col
    for i in inhibitor_list:
        unique_fixed_col_i = df.loc[df['compound_name'] == i, fixed_col].unique()
        for t1, t2 in itertools.combinations(unique_fixed_col_i, 2):
            df_i_t1 = df.loc[(df['compound_name'] == i) & (df[fixed_col] == t1)]
            df_i_t2 = df.loc[(df['compound_name'] == i) & (df[fixed_col] == t2)]
            if not (df_i_t1.empty or df_i_t2.empty):
                pairs_dict[(cell_name, i, i, t1, t2)] = pd.concat([df_i_t1, df_i_t2])

    return pairs_dict


def analyse_pairs(pairs_dict, p_value=0.05):
    keys_to_remove = []
    for key, df in pairs_dict.items():
        col_names = df.columns.tolist()
        time_col_idx = col_names.index('time')
        analysis_cols = [c for c in col_names[time_col_idx + 1:]]

        dfs_to_concat = []
        for col in analysis_cols:
            df_first = df.loc[df['compound_name'] == key[1], col]
            df_second = df.loc[df['compound_name'] == key[2], col]

            first_avg = df_first.mean()
            second_avg = df_second.mean()

            if abs(first_avg - second_avg) > p_value:
                dfs_to_concat.append(df[[col]])

        if len(dfs_to_concat) == 0:
            keys_to_remove.append(key)
        else:
            new_df = pd.concat(dfs_to_concat, axis=1)
            pairs_dict[key] = pd.concat(
                [df[['barcode', 'cell_line_name', 'compound_name', '2D_3D', 'dosage', 'time']], new_df], axis=1)

    for key in keys_to_remove:
        pairs_dict.pop(key)


def create_pairs_df(pairs_dict):
    comp_list = []
    i = 0
    for key, df in pairs_dict.items():
        cols = df.columns.tolist()
        df.iloc[0] = cols
        comp_list.append(df)
        if i < len(pairs_dict) - 1:
            comp_list.append(pd.DataFrame(np.nan, index=['-'], columns=pairs_dict[key].columns))
        i += 1

    pairs_df = pd.concat(comp_list, sort=False)
    # pd.set_option("display.max_rows", None)  # Display all rows
    # pd.set_option("display.max_columns", None)  # Display all columns
    # print(pairs_df)

    return pairs_df


def create_csv(df, name='default_name', path=None):
    """
        This method inserts the DataFrame into a csv file.

        Arguments:
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

        Arguments:
            df (pandas.DataFrame): The DataFrame to insert into a new sheet.
            path (str): The file path.
            sheet_name (str): The name of the new sheet.
    """
    if not df.empty:
        with pd.ExcelWriter(path, mode='a') as writer:
            df.to_excel(writer, sheet_name=sheet_name)
        print(f"The sheet '{sheet_name}' created successfully")
    else:
        print(f"The sheet '{sheet_name}' was not created because the DataFrame is empty")
