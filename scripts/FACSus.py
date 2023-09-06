import matplotlib.pyplot as plt
import seaborn as sns

### filters ###

# Exclude single panel sum outliers
def sanity_check_filter(df, cols, newname, lower_bound = 95, upper_bound = 105, plot = False, save = False, save_path = None):
    df[newname] = df[cols].sum(axis=1)
    n = df.shape[0]
    if plot:
        plot_hist(df, newname, 'Percentage', save = save, save_path = save_path)
    # Filtering
    df = df[(df[newname] > lower_bound) & (df[newname] < upper_bound)]
    df = df.drop(columns = [newname])
#     print('Number of samples removed for filter', newname, ':', n-df.shape[0])
    return df
# Exclude cross panel std dev outliers 
def sanity_check_filter2(df, cols, newname, plot = False, save = False, save_path = None):
    n = df.shape[0]
    df[newname] = df[cols].std(axis=1)
    df['normalised']=(df[newname]-df[newname].mean())/df[newname].std()
    if plot:
        plot_hist(df, 'normalised', 'Std Dev (Normalised)', save = save, save_path = save_path)
    # Filtering
    df = df[(df['normalised'] > -3) & (df['normalised'] < 3)]
    df[newname] = df[cols].mean(axis=1)
#     df = df.drop(columns = cols)
    df = df.drop(columns = ['normalised'])
#     print('Number of samples removed for filter', newname, ':', n-df.shape[0])
    return df
# Exclude cross panel abs diff outliers
def sanity_check_filter3(df, colname, newname, plot = False, save = False, save_path = None):
    n = df.shape[0]
    df[newname] = (df[colname + '_panel1'] - df[colname + '_panel2']).abs()
    df['normalised']=(df[newname]-df[newname].mean())/df[newname].std()
    if plot:
        plot_hist(df, 'normalised', 'Abs Diff (Normalised)', save = save, save_path = save_path)
    # Filtering
    df = df[(df['normalised'] > -3) & (df['normalised'] < 3)]
    df[newname] = (df[colname + '_panel1'] + df[colname + '_panel2'])/2
    df = df.drop(columns = ['normalised'])
    df = df.drop(columns = [colname + '_panel1', colname + '_panel2'])
#     print('Number of samples removed for filter', newname, ':', n-df.shape[0])
    return df
# Exclude single panel single marker outliers
def sanity_check_filter4(df, col, lower_bound = 3, upper_bound = -3, plot = False, save = False, save_path = None):
    n = df.shape[0]
    df['normalised']=(df[col]-df[col].mean())/df[col].std()
    if plot:
        plot_hist(df, col, 'Percentage', save = save, save_path = save_path)
    # Filtering
    df = df[(df['normalised'] > -3) & (df['normalised'] < 3)]
    df = df.drop(columns = ['normalised'])
#     print('Number of samples removed for filter', col, ':', n-df.shape[0])
    return df

### end of filters ###
### plot ###

def plot_box_distribution(df, title):
    ax = sns.boxplot(data=df, showfliers = False)
    # ax = sns.swarmplot(data=boxplot_panel1, color=".25")
    plt.title(title)
    plt.grid()
    plt.xticks(rotation=90)
    plt.ylabel('Percentage (%)')
    plt.show()
def plot_hist(df, column, xlabel, bins = 20, label = True, save = False, save_path = None):
    counts, edges, bars = plt.hist(df[column], alpha=0.75, bins = bins, histtype='bar', ec='black')
    if label:
        plt.bar_label(bars)
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title(column)
    if save:
        plt.savefig(save_path + '.png', bbox_inches = "tight", dpi=300)
    plt.show()
def plot_box(df, column):
    counts, edges, bars = plt.hist(df[column], alpha=0.75, histtype='bar', ec='black')
    plt.bar_label(bars)
    plt.xlabel('Percentage')
    plt.ylabel('Counts')
    plt.title(column)

### end of plots ###
### auxiliaries ###

def format_float(f, dec = 4):
    return f'{f:.{dec}f}'
    
### end of auxiliaries ###

### phe ###

def get_phenotype(df, col, marker_name = None, convert_decimal = True, trailing_str = 'BOSON', save_file = False):
    df['Patient ID'] = trailing_str + df['Patient ID']
    df['IID'] = df['Patient ID']
    if convert_decimal:
        df[col] = df[col]/100
    phe = df.loc[:, ['Patient ID', 'IID', col]]
    phe[col] = phe[col].apply(format_float)
    if save_file:
        if marker_name is None:
            marker_name = col
        phe.to_csv('boson_vcf/phenotypes/' + marker_name + '.txt', header = None, index = None, sep = ' ')
    return None