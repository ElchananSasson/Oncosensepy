import matplotlib.pyplot as plt
import pandas as pd
import os
import exceptions as e


def is_valid_L(df):
    """
        This method checks whether the DataFrame is in the appropriate format.

        Arguments:
            df (pandas.DataFrame): The DataFrame to be checked.

        Returns:
            bool: True if valid, throw an exception otherwise.
    """
    if df.columns[0] != 'barcode':
        raise e.InvalidDataSetException("column 0 should be 'barcode'")
    if df.columns[1] != 'cell_line_name':
        raise e.InvalidDataSetException("column 1 should be 'cell_line_name'")
    if df.columns[2] != 'compound_name':
        raise e.InvalidDataSetException("column 2 should be 'compound_name'")
    if df.columns[3] != 'D2_D3':
        raise e.InvalidDataSetException("column 3 should be 'D2_D3'")
    if df.columns[4] != 'dosage':
        raise e.InvalidDataSetException("column 3 should be 'dosage'")
    if df.columns[5] != 'time':
        raise e.InvalidDataSetException("column 4 should be 'time'")
    return True


def is_valid_path(path, directory=True):
    """
        This method checks whether the DataFrame is in the appropriate format.

        Arguments:
            path (str): The path to be checked.
            directory (bool): Default ia True. If the path is directory - fill True, otherwise - False.

        Returns:
            bool: True if valid, throw an exception otherwise.
    """
    if not os.path.exists(path):
        raise e.InvalidPathException(f"The path '{path}' didn't exist")
    if directory:
        if not os.path.isdir(path):
            raise e.InvalidDirectoryPathException(f"The path '{path}' should be to directory")
    return True


def is_valid_G(df):
    """
        This method checks whether the DataFrame is in the appropriate format.

        Arguments:
            df (pandas.DataFrame): The DataFrame to be checked.

        Returns:
            bool: True if valid, throw an exception otherwise.
    """
    if df.columns[0] != 'UID':
        raise e.InvalidDataSetException("column 0 should be 'UID'")

    return True


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
    is_valid_L(df)
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

        Arguments:
            df (pandas.DataFrame): The DataFrame to filter.
            col (str): The column according to which the filtering will be performed.
            filter_list (List): A list of values that we would like to appear in the selected column.

        Returns:
            pandas.DataFrame: The DataFrame after filtering.
    """
    is_valid_L(df)
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
    is_valid_path(path)
    plt.title(title)
    plt.figure(figsize=(25, 15))
    plt.scatter(uid, values)
    plt.xticks(uid, [f'{name} ({i})' for i, name in enumerate(uid)], rotation=90, fontsize=6)
    plt.ylabel('Effect')

    for i, txt in enumerate(range(len(uid))):
        plt.text(uid[i], values[i] + 0.002, str(i), fontsize=5)

    plt.savefig(path + '/' + f'{title}.jpg', dpi=300)
    plt.close()


def sort_plot_G_values(g_df, cols, path):
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
    is_valid_path(path)
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
