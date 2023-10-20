import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import numpy as np
# from scipy import interpolate

# Data Read
def read_INPE(fileName, alphaINPE):

    dfINPE = pd.read_csv('Pontos\\'+fileName+'.txt', skiprows=12, names=['ABS_Ant', 'PHA_Ant', 'ABS_Ref', 'PHA_Ref', 'E'])
    dfINPE = dfINPE.iloc[:-1, :-1]
    dfINPE['ABS_Ant'] = pd.to_numeric(dfINPE['ABS_Ant'], errors='coerce')
    # dfINPE['POSr'] = np.radians(range(360+alphaINPE, 0+alphaINPE,-1))
    config = pd.read_csv('Pontos\\'+fileName+'.txt', skiprows=2, nrows=3, names = ['#','npts', 'start', 'stop', 'roll', 'direction','0']).iloc[0]

    stepROLL = -1 if config['direction']== 'REVERSE'else 1

    dfINPE['POSr'] = np.radians(range(int(config['start'])+alphaINPE,int(config['stop'])+alphaINPE+stepROLL, stepROLL))

    return dfINPE

def read_CST(fileName, alphaCST):
    
    dfCST = pd.read_csv('CST\\'+fileName+'.txt', sep='\s+', skiprows=2, names= ['Theta [deg.]','Phi [deg.]','Abs(Grlz)[dBi]','Abs(Theta)[dBi]',  'Phase(Theta)[deg.]', 'Abs(Phi)[dBi]', 'Phase(Phi)[deg.]', 'Ax.Ratio[dB]'])
    dfCST['POSr'] = np.radians(range(0+alphaCST, 360+alphaCST, 1))

    return dfCST

# Normalize function
def normalize_to_Ref (df1,col_to_normalize, df2,col_refence):

    max_df2 = df2[col_refence].max()
    df1[col_to_normalize] = df1[col_to_normalize] +(max_df2 - df1[col_to_normalize].max())

    return df1

def normalize_to_0 (df1,col1, df2,col2):

    df1[col1] = df1[col1] +abs(0 - df1[col1].max()) if df1[col1].max() < 0 else df1[col1] -abs(0 - df1[col1].max()) 
    df2[col2] = df2[col2] +abs(0 - df2[col2].max()) if df2[col2].max() < 0 else df2[col2] -abs(0 - df2[col2].max())

    return df1, df2

# Correlation
def correlation(df1, col1, df2, col2):
    d= { 'INPE':df1[col1],'CST':df2[col2]}

    df = pd.DataFrame(d)
    corr_pearson = df.corr(method='pearson')
    corr_kendall = df.corr(method='kendall')
    corr_spearman = df.corr(method='spearman')

    return corr_pearson, corr_kendall, corr_spearman

#PLOT
def plot (dfINPE, dfCST, fileName, cutAngle):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='polar')
    ax.plot(dfINPE['POSr'], dfINPE['ABS_Ant'], label='INPE', color='red')
    ax.plot(dfCST['POSr'], dfCST['Abs('+cutAngle+')[dBi]'], '--', label='CST', color='blue',)

    plt.yticks(range(round(max(dfCST['Abs('+cutAngle+')[dBi]'])-30), round(max(dfCST['Abs('+cutAngle+')[dBi]'])), 3))
    # ax.set(ylim=(max(dfCST['Abs(Theta)[dBi]'])-20, max(dfCST['Abs(Theta)[dBi]'])))
    ax.set_title('Radiation Diagram - '+fileName)
    ax.legend()
    
    plt.savefig('results/'+fileName+'.png', dpi=300, bbox_inches='tight')
    plt.savefig('results/'+fileName+'.svg', format='svg', bbox_inches='tight')
    
    plt.show()

    return 

if __name__ == "__main__":
    #PARAMETERS
    antenna ='2A'
    tx = 'V'
    diag = 'V'
    freq = 'F09'
    fileName = antenna+tx+diag+freq
    cutAngle = 'Theta' if diag == 'V' else 'Phi' # H

    #ROTATION
    alphaINPE = 0
    alphaCST = 90

    # INPE DATA
    dfINPE = read_INPE (fileName, alphaINPE)
    
    # CST DATA
    dfCST = read_CST(fileName, alphaCST)

    # dfINPE = normalize_to_Ref(dfINPE,'ABS', dfCST,'Abs('+cutAngle+')[dBi]')
    dfINPE, dfCST = normalize_to_0(dfINPE,'ABS_Ant', dfCST,'Abs('+cutAngle+')[dBi]')

    corr_pearson, corr_kendall, corr_spearman = correlation(dfINPE,'ABS_Ant', dfCST,'Abs('+cutAngle+')[dBi]')
    print (corr_pearson, corr_kendall, corr_spearman)

    # PLOT
    plot(dfINPE, dfCST, fileName, cutAngle)