import validation as valid
import matplotlib.pyplot as plt


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
