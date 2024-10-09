import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({'font.size': 10})
plt.rcParams.update({'figure.autolayout': True})
# import healpy as hp
# import pymaster as nmt

comparison_group = 2
JOB = '/mnt/d/data_large/unwise_sz/'
PATHRESULTS = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/'

##
# original txt writing
# np.savetxt(OUTPATH + f'namaster_comparison{comparison_group}.txt',np.column_stack(cols),header=' '.join(names))
##
# Load data and header, generate table structure

data = pd.read_csv(PATHRESULTS + f'namaster_comparison{comparison_group}.txt', sep=' ')
colum_list = data.columns.tolist()
values = data.values

table = {}
for i in range(len(colum_list)):
    table[colum_list[i]] = values[:,i]

if comparison_group == 0:
    #CIB deprojection
    ymap_name_list = ['no_deprojection',
                        'CIB+CMB_T=10.17beta=1.7',
                        'CIB+CMB_T=24beta=1.0',
                        'CIB+CMB_T=24beta=1.4',
                        'CIB+CMB_T=10.14beta=1.4',
                        'CIB+CMB_T=10.14beta=1.6']

    ymap_path_list = [
                    JOB+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/deproject_CMB_CIB_default_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.0_T24_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.4_T24_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.4_T10.14_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.6_T10.14_standard_full.fits']
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
                    JOB+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/deproject_CMB5_CIB_CIBdbeta_default_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.0_T24_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.4_T24_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.4_T10.14_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.6_T10.14_standard_full.fits']
elif comparison_group == 2:
    ymap_name_list = [
                    'no_deprojection',
                    'CMB+CIB+CIBdbetadT_T=10.17beta=1.7',
                    'CMB+CIB+CIBdbetadT_T=24beta=1.0',
                    'CMB+CIB+CIBdbetadT_T=24beta=1.4',
                    'CMB+CIB+CIBdbetadT_T=10.14beta=1.4',
                    'CMB+CIB+CIBdbetadT_T=10.14beta=1.6']

    ymap_path_list = [
                    JOB+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/deproject_CMB5_CIB_CIBdbeta_CIBdT_default_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.0_T24_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.4_T24_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.4_T10.14_standard_full.fits',
                    JOB+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.6_T10.14_standard_full.fits']


# IMGPATH = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/images/'

# # print first few lines
# # print(data[:5])

# omit = 4
# truncate = 50
# fig,ax = plt.subplots(1,2,figsize=(12,6))
# ell = table['ell'][omit:truncate]
# pll = ell*(ell+1)/(2*np.pi)

# for i in range(len(ymap_name_list)):
#     ax[0].plot(ell, pll*table[ymap_name_list[i]+'_yy'][omit:truncate],label=ymap_name_list[i])
#     ax[1].plot(ell, pll*table[ymap_name_list[i]+'_gy'][omit:truncate],label=ymap_name_list[i])

# ax[0].legend()
# ax[0].set_ylabel(r'$\ell(\ell+1)C_{\ell}^{yy}/2\pi$')
# ax[0].set_xlabel(r'$\ell$')
# ax[0].set_title('Cl_yy')
# ax[1].legend()
# ax[1].set_ylabel(r'$\ell(\ell+1)C_{\ell}^{gy}/2\pi$')
# ax[1].set_xlabel(r'$\ell$')
# ax[1].set_title('Cl_gy')
# fig.savefig(IMGPATH+f'CIBdepro_comparison{comparison_group}.png',dpi=400)
