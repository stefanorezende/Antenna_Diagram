import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import numpy as np
# from scipy.stats import pearsonr

# Data Read
def read_INPE(fileName, alphaINPE):

    dfINPE = pd.read_csv('Pontos\\'+fileName+'.txt', skiprows=12, names=['ABS_Ant', 'PHA_Ant', 'ABS_Ref', 'PHA_Ref', 'E'])
    dfINPE = dfINPE.iloc[:-1, :-1]
    dfINPE['ABS_Ant'] = pd.to_numeric(dfINPE['ABS_Ant'], errors='coerce')
    # dfINPE['POSr'] = np.radians(range(360+alphaINPE, 0+alphaINPE,-1))
    config = pd.read_csv('Pontos\\'+fileName+'.txt', skiprows=2, nrows=3, names = ['#','npts', 'start', 'stop', 'roll', 'direction','0']).iloc[0]

    stepROLL = -1 if config['direction']== 'REVERSE'else 1

    # dfINPE['POSr'] = np.radians(range(int(config['start'])+alphaINPE,int(config['stop'])+alphaINPE+stepROLL, stepROLL))
    dfINPE['POS'] = range(int(config['start'])+alphaINPE,int(config['stop'])+alphaINPE+stepROLL, stepROLL)

    dfINPE.loc[dfINPE['POS'] < 0, 'POS'] += 360
    dfINPE = dfINPE.sort_values(by= 'POS', ascending=True)

    df180 = dfINPE.head(180)

    dfINPE = dfINPE.iloc[180:]

    dfINPE = pd.concat([dfINPE, df180], ignore_index=True)

    dfINPE ['POSr'] = np.radians(dfINPE['POS'])

    return dfINPE

def read_CST(fileName, alphaCST):
    
    dfCST = pd.read_csv('CST\\'+fileName+'.txt', sep='\s+', skiprows=2, names= ['Theta [deg.]','Phi [deg.]','Abs(Grlz)[dBi]','Abs(Theta)[dBi]',  'Phase(Theta)[deg.]', 'Abs(Phi)[dBi]', 'Phase(Phi)[deg.]', 'Ax.Ratio[dB]'])
    # dfCST['POSr'] = np.radians(range(0+alphaCST, 360+alphaCST, 1))

    # df90 = dfCST.head(90)

    # dfCST = dfCST.iloc[90:]

    # dfCST = pd.concat([dfCST, df90], ignore_index=True)
    dfCST['POSr'] = np.radians(range(0, 360))
    return dfCST

# GAIN CALC
def gain(diag, freq, df1,col1, df2,col2):
    sAnt = df1[col1].max()

    sHorndf = pd.read_csv('AGP\\P'+diag+freq+'.txt', skiprows=12, names=['ABS_Ant', 'PHA_Ant', 'ABS_Ref', 'PHA_Ref', 'E'])
    sHorndf = sHorndf.iloc[:-1, :-1]
    sHorndf['ABS_Ant'] = pd.to_numeric(sHorndf['ABS_Ant'], errors='coerce')
    sHorn = sHorndf['ABS_Ant'].max()

    freqInt= int(freq[-2:])*10**9

    gHorndf = pd.read_excel('Ganho.xlsx', sheet_name='AGP')
    gHorn = float(gHorndf.loc[gHorndf['ID']== freq, 'Ganho [dBi]'])

    prGdf = pd.read_csv('PR\\Cabo_Gore.prn',skiprows=1)
    prGdf = prGdf.iloc[:, :-1]

    prHdf = pd.read_csv('PR\\Cabo_Huber.prn',skiprows=1)
    prHdf = prHdf.iloc[:, :-1]

    fprG = float(prGdf.loc[prGdf['Frequency (Hz)']== freqInt, 'dB'])
    fprH = float(prHdf.loc[prHdf['Frequency (Hz)']== freqInt, 'dB'])

    gINPE = gHorn - ((sHorn - sAnt) - (fprH - fprG))+6.2

    gCST = df2[col2].max()

    return gCST, gINPE

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
    # corr_pearson = pearsonr(df['INPE'],df['CST'])
    corr_pearson = df['CST'].corr(df['INPE'], method='pearson')
    corr_kendall = df['CST'].corr(df['INPE'], method='kendall')
    corr_spearman = df['CST'].corr(df['INPE'], method='spearman')

    return corr_pearson, corr_kendall, corr_spearman

#PLOT
def plot (dfINPE, dfCST, fileName, cutAngle,corr_pearson,gCST, gINPE):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='polar')
    ax.plot(dfINPE['POSr'], dfINPE['ABS_Ant'], label='Meas',  color='red')
    ax.plot(dfCST['POSr'], dfCST['Abs('+cutAngle+')[dBi]'], '--', label='Simul', color='blue',)  
    
    ax.set_theta_offset(np.pi / 2)
    # plt.yticks(range(round(max(dfCST['Abs('+cutAngle+')[dBi]'])-30), round(max(dfCST['Abs('+cutAngle+')[dBi]'])), 10))
    # ax.set(ylim=(max(dfCST['Abs(Theta)[dBi]'])-30, max(dfCST['Abs(Theta)[dBi]'])))
    ax.set (ylim = (-40, 0))
    ax.set_title('Radiation Diagram')
    ax.legend()
    text = 'G_simul: '+str(round(gCST,1))+' dBi'+ '\nG_meas: '+ str(round(gINPE,1)) + ' dBi'
    # plt.text( 9.5,9.5,text , fontsize=10, bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round'))
    plt.text(0.8, -0.05, text, transform=ax.transAxes, fontsize=10, verticalalignment='bottom', bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round'))
    plt.text(0.05, 1, 'Freq: 11 GHz', transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round'))
    plt.savefig('results\\'+fileName+'.png', dpi=300, bbox_inches='tight')
    plt.savefig('results\\'+fileName+'.svg', format='svg', bbox_inches='tight')
    
    plt.show()

    return 

def plotlin(dfINPE, dfCST, fileName, cutAngle):
    degree = list(np.deg2rad(range(0,360)))
    
    figure(figsize=(15, 12), dpi=300)
    fig, ax = plt.subplots( layout = 'constrained')
    ax.plot(degree, dfINPE['ABS_Ant'], label='INPE', color='red')
    ax.plot(degree, dfCST['Abs('+cutAngle+')[dBi]'], '--', label='CST', color='blue',)

    ax.set_xlabel('Frequency')  # Add an x-label to the axes.
    ax.set_ylabel('S(1,1) [dB]')  # Add a y-label to the axes.
    ax.set_title('S11_Analysis')  # Add a title to the axes.

    ax.grid()
    ax.legend()
    plt.show()

    return


if __name__ == "__main__":
    #PARAMETERS
    antenna ='1A'
    tx = 'H'
    diag = 'H'
    freq = 'F11'
    fileName = antenna+tx+diag+freq
    cutAngle = 'Theta' if diag == 'V' else 'Phi' # H

    #ROTATION
    alphaINPE = 50

    # alphaCST = -90 if diag == 'H' else 90
    alphaCST = 0

    # INPE DATA
    dfINPE = read_INPE (fileName, alphaINPE)
    
    # CST DATA
    dfCST = read_CST(fileName, alphaCST)

    gCST, gINPE = gain(diag, freq, dfINPE,'ABS_Ant', dfCST,'Abs('+cutAngle+')[dBi]')
    
    # dfINPE = normalize_to_Ref(dfINPE,'ABS', dfCST,'Abs('+cutAngle+')[dBi]')
    dfINPE, dfCST = normalize_to_0(dfINPE,'ABS_Ant', dfCST,'Abs('+cutAngle+')[dBi]')

    corr_pearson, corr_kendall, corr_spearman = correlation(dfINPE,'ABS_Ant', dfCST,'Abs('+cutAngle+')[dBi]')
    
    with open('results\\'+fileName+'_corr.txt', 'w') as arquivo:
        arquivo.write(f'Pearson: {corr_pearson}\nKendall: {corr_kendall}\nSpearman: {corr_spearman}')
        arquivo.write(f'\nGain CST: {gCST}\nGain INPE: {round(gINPE,2)}')

    # PLOT
    plot(dfINPE, dfCST, fileName, cutAngle, corr_pearson, gCST, gINPE)
    # plotlin(dfINPE, dfCST, fileName, cutAngle)