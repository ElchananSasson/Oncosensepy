import pandas as pd
import data_science as ds

if __name__ == '__main__':
    data_name = 'Table1_myData66'
    data_set_path = r'Data/' + data_name + '.xlsx'

    l_df = pd.read_excel(data_set_path, sheet_name='L').fillna(0)
    g_df = pd.read_excel(data_set_path, sheet_name='G').fillna(0)
    err_limit_lambda = pd.read_excel(data_set_path, sheet_name='ErrorLimitLambda').columns.values[0]

    important_l = ds.important_L(l_df, err_limit_lambda, 2)
    # ds.create_new_sheet(important_l, data_set_path, 'important_L')
    #
    # filter_dosage = ds.filter_by_col(important_l, 'dosage', ['0.00001nm', '1nm', '40nm', '1uM'])
    # ds.create_new_sheet(filter_dosage, data_set_path, 'filter_by_dosage')
    #
    # filter_time = ds.filter_by_col(filter_dosage, 'time', ['0hr', '24hr'])
    # ds.create_new_sheet(filter_time, data_set_path, 'filter_by_dosage_time')
    #
    important_g = ds.sort_plot_G_values(g_df, important_l.columns[6:], r'G_plot')
    ds.create_new_sheet(important_g, data_set_path, 'Sorted_G')

    # comp = ds.compare_pairs(important_l, 'MCF7', ['CONTROL', 'DMSO', 'GSK1838705A', 'LAPATINIB'])
    # ds.create_new_sheet(comp, data_set_path, 'MCF7_compare')
