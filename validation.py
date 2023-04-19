import exceptions as e
import os


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
    if df.columns[3] != '2D_3D':
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
