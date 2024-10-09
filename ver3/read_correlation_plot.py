import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTDIR = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/'
IMGDIR = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/images'
sample = 'nobeam'
result_gy_path = RESULTDIR + sample + '/namaster_comparison_gy.csv'
result_yy_path = RESULTDIR + sample + '/namaster_comparison_yy.csv'
result_gg_path = RESULTDIR + sample + '/namaster_comparison_gg.csv'

df_gy = pd.read_csv(result_gy_path)
df_yy = pd.read_csv(result_yy_path)
df_gg = pd.read_csv(result_gg_path)

ell = df_gy['ell'].values
gg = df_gg['cl_gg'].values

header_gy = df_gy.columns
header_yy = df_yy.columns

comparison_group = 0

if comparison_group == 0:
    #CIB deprojection
    ymap_name_list = [
                    'no_deprojection',
                    'CIB+CMB_T=10.17beta=1.7',
                    'CIB+CMB_T=24beta=1.0',
                    'CIB+CMB_T=24beta=1.4',
                    'CIB+CMB_T=10.14beta=1.4',
                    'CIB+CMB_T=10.14beta=1.6']

elif comparison_group == 1:
    #CIB + CIBdbeta deprojection
    ymap_name_list = [
                    'no_deprojection',
                    'CMB+CIB+CIBdbeta_T=10.17beta=1.7',
                    'CMB+CIB+CIBdbeta_T=24beta=1.0',
                    'CMB+CIB+CIBdbeta_T=24beta=1.4',
                    'CMB+CIB+CIBdbeta_T=10.14beta=1.4',
                    'CMB+CIB+CIBdbeta_T=10.14beta=1.6']

elif comparison_group == 2:
    ymap_name_list = [
                    'no_deprojection',
                    'CMB+CIB+CIBdbetadT_T=10.17beta=1.7',
                    'CMB+CIB+CIBdbetadT_T=24beta=1.0',
                    'CMB+CIB+CIBdbetadT_T=24beta=1.4',
                    'CMB+CIB+CIBdbetadT_T=10.14beta=1.4',
                    'CMB+CIB+CIBdbetadT_T=10.14beta=1.6']

ploystyles = ['-', '--', ':', ':', ':', ':']
BIN = 50
f_sky = 0.572050134340922


plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'figure.autolayout': True})

fig,ax = plt.subplots(1,2,figsize=(12,6))
omit = 4

# set thicker border for figure
for axis in ['top','bottom','left','right']:
    ax[0].spines[axis].set_linewidth(1.5)
    ax[1].spines[axis].set_linewidth(1.5)

for model_name in ymap_name_list:
    model_name_gy = model_name + '_gy'
    model_name_yy = model_name + '_yy'
    gy = df_gy[model_name_gy].values
    yy = df_yy[model_name_yy].values
    
    sigma_gy = np.sqrt((gg*yy+gy**2)/((2*ell+1)*BIN*f_sky))
    sigma_yy = np.sqrt((2*yy**2)/((2*ell+1)*BIN*f_sky))
    
    pll = ell*(ell+1)/(2*np.pi)
    
    Cl_gy = gy*pll
    Cl_yy = yy*pll
    Sigma_gy = sigma_gy*pll
    Sigma_yy = sigma_yy*pll
    
    ax[0].errorbar(ell[omit:], Cl_gy[omit:], yerr=Sigma_gy[omit:], label=model_name, linestyle=ploystyles[ymap_name_list.index(model_name)])
    ax[1].errorbar(ell[omit:], Cl_yy[omit:], yerr=Sigma_yy[omit:], label=model_name, linestyle=ploystyles[ymap_name_list.index(model_name)])

ax[1].set_title('Comparison of $C_{\ell}^{yy}$')
ax[0].set_xticks([10, 100, 1000])
ax[0].set_xscale('log')
ax[0].set_xlabel(r'$\ell$')
ax[0].set_ylabel(r'$\ell(\ell+1)C_{\ell}^{gy}/2\pi$')
ax[0].legend(fontsize = 'small')

ax[1].set_title('Comparison of $C_{\ell}^{yy}$')
ax[1].set_xticks([10, 100, 1000])
ax[1].set_xscale('log')
ax[1].set_xlabel(r'$\ell$')
ax[1].set_ylabel(r'$\ell(\ell+1)C_{\ell}^{yy}/2\pi$')
ax[1].legend(fontsize = 'small')

fig.savefig(IMGDIR + f'/namaster_comparison{comparison_group}_sample{sample}.png')
    

