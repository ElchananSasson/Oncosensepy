import pandas as pd
import os


class InvalidDataSetException(Exception):
    def __init__(self, message):
        super().__init__(message)


class LambdaDataSet:
    def __init__(self):
        self.columns = None

    @staticmethod
    def is_valid_L_data_set(df):
        """
            This method checks whether the DataFrame is in the appropriate format.

            Arguments:
                df (pandas.DataFrame): The DataFrame to be checked.

            Returns:
                bool: True if valid, throw an exception otherwise.
        """
        if df.columns[0] != "barcode":
            raise InvalidDataSetException("column 0 should be 'barcode'")
        if df.columns[1] != "cell_line_name":
            InvalidDataSetException("column 1 should be 'cell_line_name'")
        if df.columns[2] != "compound_name":
            InvalidDataSetException("column 2 should be 'compound_name'")
        if df.columns[3] != "dosage":
            InvalidDataSetException("column 3 should be 'dosage'")
        if df.columns[4] != "time":
            InvalidDataSetException("column 4 should be 'time'")
        return True

    def important_L(self, df, err_limit, threshold):
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
        if self.is_valid_L_data_set(df):
            new_df = df.loc[:, :4].copy()
            for i in range(1, len(df.columns) - 4):
                count = 0
                for cell in df[i]:
                    if abs(cell) > err_limit:
                        count += 1
                if count >= threshold:
                    new_df[i] = df[i].copy()
            return new_df

    def filter_by_col(self, df, col, filter_list):
        """
            This function filters data by a certain column and by a list of values it receives.

            Arguments:
                df (pandas.DataFrame): The DataFrame to filter.
                col (str): The column according to which the filtering will be performed.
                filter_list (List): A list of values that we would like to appear in the selected column.

            Returns:
                pandas.DataFrame: The DataFrame after filtering.
        """
        if self.is_valid_L_data_set(df):
            filter_df = df.loc[(df[col].isin(filter_list))]
            if len(filter_df) == 0:
                print(f"There is no data to show by '{col}' filtering")
            return filter_df


def sort_G(g_df, df):
    cols = df.columns[5:]
    important_g = pd.DataFrame(columns=pd.MultiIndex.from_product([cols, ['UID', 'Effect']]), index=False)
    for col, sub in important_g.columns:
        if sub == 'UID':
            important_g[(col, sub)] = g_df['UID'].copy()
        else:
            important_g[(col, sub)] = g_df[col]
        print(important_g)
        # return important_g
        # important_g[col] = g_df[col[]]
    # temp = pd.DataFrame()
    # temp = pd.DataFrame(columns=pd.MultiIndex.from_product([[1], ['UID', 'Effect']]))
    # for col in important_g.columns:
    #     if col[1] == 'UID':
    #         temp[(1, 'UID')] = g_df['UID']
    #     if col[1] == 'Effect':
    #         temp[(1, 'Effect')] = g_df[col[0]]
    #         print(temp)
    #         temp = temp.sort_values(by=temp[(1, 'UID')].apply(lambda x: x[(1, 'Effect')]), inplace=True)
    #         print(temp)
    #         important_g[col] = temp
    #         print(important_g[col])
    #         temp = pd.DataFrame(columns=pd.MultiIndex.from_product([[1], ['UID', 'Effect']]))
            # important_g[col] = g_df[col[0]]
            # important_g[col[0]] = important_g[col[0]].sort_values(by=col[1])

            # print(temp[(col[0], 'UID')], temp[(col[0], 'Effect')])
            # important_g[col[0]] = temp
            # print(important_g[col[0]])
            # important_g[(col[0], 'Effect')] = temp[(col[0], 'Effect')]
            # print(important_g[col[0]])
            # temp = pd.DataFrame()

    return important_g


def create_csv(df, name='default_name', path='/'):
    """
        This method inserts the DataFrame into a csv file.

        Arguments:
            df (pandas.DataFrame): The DataFrame to insert into a new csv.
            name (str): File name.
            path (str): The file path.
    """
    if df.empty:
        return
    # Create the folder if it does not exist
    if path == '/':
        folder = 'csv'
        if not os.path.exists(folder):
            os.makedirs(folder)
        path = folder + '/' + name + '.csv'  # Save the data frame to the CSV file
        df.to_csv(path, index=False)
    else:
        df.to_csv(path + '/' + name + '.csv', index=False)
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
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        print(f"The sheet '{sheet_name}' created successfully")
    else:
        print(f"The sheet '{sheet_name}' was not created because the DataFrame is empty")
