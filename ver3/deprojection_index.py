
DAT = '/mnt/d/data_large/unwise_sz/'
import os

def get_ymap_index_planck(comparison_group):
    if comparison_group == 0:
        #CIB deprojection
        ymap_name_list = ['no_deprojection',
                            'CIB+CMB_T=10.17beta=1.7',
                            'CIB+CMB_T=24beta=1.0',
                            'CIB+CMB_T=24beta=1.4',
                            'CIB+CMB_T=10.14beta=1.4',
                            'CIB+CMB_T=10.14beta=1.6']

        ymap_path_list = [
                        DAT+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/deproject_CMB_CIB_default_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.0_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.4_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.4_T10.14_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.6_T10.14_standard_full.fits']
    elif comparison_group == 1:
        #CIB + CIBdbeta deprojection
        ymap_name_list = [
                        'no_deprojection',
                        'CMB+CIB+CIBdbeta_T=10.17beta=1.7',
                        'CMB+CIB+CIBdbeta_T=24beta=1.0',
                        'CMB+CIB+CIBdbeta_T=24beta=1.4',
                        'CMB+CIB+CIBdbeta_T=10.14beta=1.4',
                        'CMB+CIB+CIBdbeta_T=10.14beta=1.6']
        ymap_path_list = [
                        DAT+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/deproject_CMB5_CIB_CIBdbeta_default_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.0_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.4_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.4_T10.14_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.6_T10.14_standard_full.fits']
    elif comparison_group == 2:
        ymap_name_list = [
                        'no_deprojection',
                        'CMB+CIB+CIBdbetadT_T=10.17beta=1.7',
                        'CMB+CIB+CIBdbetadT_T=24beta=1.0',
                        'CMB+CIB+CIBdbetadT_T=24beta=1.4',
                        'CMB+CIB+CIBdbetadT_T=10.14beta=1.4',
                        'CMB+CIB+CIBdbetadT_T=10.14beta=1.6']

        ymap_path_list = [
                        DAT+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/deproject_CMB5_CIB_CIBdbeta_CIBdT_default_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.0_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.4_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.4_T10.14_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.6_T10.14_standard_full.fits']
    
    return ymap_name_list, ymap_path_list


class deprojection_index_act:
    def __init__(self,directory:str,filename:str):
        self.dir = directory
        self.filename = filename
        self.deprotype = filename.split('_')[-3]
        self.nu = float(filename.split('_')[-2])
        self.T = float(filename.split('_')[-1][:-5])
        self.path = directory+filename

def get_ymap_index_act():
    filename_cib = os.listdir(DAT+'ACT/deprojections/cib/')
    filename_cibdbeta = os.listdir(DAT+'ACT/deprojections/cib_cibdbeta/')
    
    codex = []
    for i in range(len(filename_cib)):
        index =  deprojection_index_act(DAT+'ACT/deprojections/cib/',filename_cib[i])
        codex.append(index)
    for i in range(len(filename_cibdbeta)):
        index =  deprojection_index_act(DAT+'ACT/deprojections/cib_cibdbeta/',filename_cibdbeta[i])
        codex.append(index)
    return codex

if __name__ == '__main__':
    codex = get_ymap_index_act()
    
    for i in range(len(codex)):
        print(codex[i].deprotype,codex[i].nu,codex[i].T)