import pandas as pd
import oncosensepy as osp

if __name__ == '__main__':
    data_name = 'Table1_myData66'
    data_set_path = r'Data/' + data_name + '.xlsx'

    l_df, g_df = osp.get_LG_data(data_set_path)
    err_limit_lambda = pd.read_excel(data_set_path, sheet_name='ErrorLimitLambda').columns.values[0]

    important_l = osp.important_L(l_df, err_limit_lambda, 2)
    # osp.create_new_sheet(important_l, data_set_path, 'important_L')

    # filter_dosage = osp.filter_by_col(important_l, 'dosage', ['0.00001nm', '1nm', '40nm', '1uM'])
    # osp.create_new_sheet(filter_dosage, data_set_path, 'filter_by_dosage')

    # filter_time = osp.filter_by_col(filter_dosage, 'time', ['0hr', '24hr'])
    # osp.create_new_sheet(filter_time, data_set_path, 'filter_by_dosage_time')

    # important_g = osp.sort_G_values(g_df, important_l.columns[6:], '', save=False)
    # osp.create_new_sheet(important_g, data_set_path, 'Sorted_G')

    # pairs_dict = osp.pairs_df_to_dict(important_l, 'MDAMB468', fixed_col='time')  # data27
    pairs_dict = osp.pairs_df_to_dict(important_l, 'MCF7', fixed_col='time')  # data66
    # pairs_dict = osp.pairs_df_to_dict(important_l, 'MCAS', fixed_col='time')  # data69

    analyzed_pairs_dict = osp.analyze_pairs(pairs_dict, 0.05, display=False, only_avg=True)
    pairs_df = osp.create_pairs_df(analyzed_pairs_dict)
    osp.create_new_sheet(pairs_df, data_set_path, 'MCF7')
    #
    # control_treatment_df = osp.analyze_control_treatment(important_l, 'MCF7')
    # osp.create_new_sheet(control_treatment_df, data_set_path, 'MCF7_control_vs_treatment')
