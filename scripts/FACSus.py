import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import scipy as sp
import csv
import io
import os
from scipy.stats import poisson
import itertools
import collections
import seaborn as sns
from sklearn import preprocessing
from PIL import Image, ImageDraw, ImageFont
from IPython.display import display
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
        plot_strip(df, col, save = save, save_path = save_path)
    # Filtering
    df = df[(df['normalised'] > -3) & (df['normalised'] < 3)]
    df = df.drop(columns = ['normalised'])
#     print('Number of samples removed for filter', col, ':', n-df.shape[0])
    return df

### end of filters ###
### plot ###

def plot_strip(df, column, xlabel = None, title = None, upper = 0.75, lower = 0.25, save = False, save_path = None):
    data = df[column]
    upper_quartile = data.quantile(upper)
    lower_quartile = data.quantile(lower)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    sns.stripplot(data=data, jitter = 0.01, ax = ax, alpha = 0.6, zorder = 1)
    ax.hlines(upper_quartile, xmin=-0.05, xmax=0.05, colors="black", linestyles="-", linewidth=2, alpha = 0)
    ax.hlines(lower_quartile, xmin=-0.05, xmax=0.05, colors="black", linestyles="-", linewidth=2, alpha = 0)
    ax.hlines(upper_quartile, xmin=-0.005, xmax=0.005, colors="black", linestyles="-", linewidth=2)
    ax.hlines(lower_quartile, xmin=-0.005, xmax=0.005, colors="black", linestyles="-", linewidth=2)
    ax.vlines(x = 0, ymin = upper_quartile, ymax = lower_quartile, colors = "black", linestyles = "-", linewidth = 2)
    ax.scatter(0, data.mean(), marker = 'x', alpha = 1, s = 40, c = 'black')
    
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.set_ylabel('Percentage (%)')
    if save:
        plt.savefig(save_path, bbox_inches = "tight", dpi=300)
def plot_hist(df, column, xlabel = None, bins = 20, label = True, save = False, save_path = None):
    counts, edges, bars = plt.hist(df[column], alpha=0.75, bins = bins, histtype='bar', ec='black')
    if label:
        plt.bar_label(bars)
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title(column)
    if save:
        plt.savefig(save_path, bbox_inches = "tight", dpi=300)
    plt.show()
def plot_box(df, column):
    counts, edges, bars = plt.hist(df[column], alpha=0.75, histtype='bar', ec='black')
    plt.bar_label(bars)
    plt.xlabel('Percentage')
    plt.ylabel('Counts')
    plt.title(column)
def combine_gwas_res(marker, indir = "/gpfs3/well/ansari/users/gjx698/BOSON_FACS/graphs/facs-host/baseline/all/", 
                     outdir = None, title_size = 200, show = False, save = False,
                     font_path = "/well/ansari/users/gjx698/conda/skylake/envs/sus/fonts/SourceCodePro-BlackIt.ttf", 
                     output_resize = True, output_resize_multiplier = 5):
    image1 = Image.open(indir + "Rect_Manhtn." + marker + ".jpg")
    image2 = Image.open(indir + "strip." + marker + ".png")
    image3 = Image.open(indir + "QQplot." + marker + ".jpg")
    
    # Get the dimensions of the images
    width1, height1 = image1.size
    width2, height2 = image2.size
    width3, height3 = image3.size
    
    if height2 >= height3:
        height3 = height2
        width3 = int(width3 * height2/height3)
        image3.resize((width3, height3))
    else:
        height2 = height3
        width2 = int(width2 * height3/height2)
        image2.resize((width2, height2))
    
    # Calculate the width and height of the combined image
    total_width = width1
    max_height = height1 + height2 + title_size
    
    # Create a new blank image with the calculated dimensions
    combined_image = Image.new("RGB", (total_width, max_height), (255, 255, 255))
    
    interval = width1 - width2 - width3
    blank = int(interval/3)
    
    # Paste the three images onto the combined image
    combined_image.paste(image1, (0, height2+title_size))
    combined_image.paste(image2, (blank, title_size + 150))
    combined_image.paste(image3, (width2 + blank*2, title_size))
    
    draw = ImageDraw.Draw(combined_image)
    
    # Define the title or caption text and font
    
    font = ImageFont.truetype(font_path, size = title_size)  # You can specify a different font and size
    
    # Define the position where you want to place the title
    title_position = (0, 10)  # Adjust the coordinates as needed
    
    # Define the text color (RGB format)
    text_color = (0, 0, 0)  # Black color
    
    draw.text(title_position, marker, fill=text_color, font=font)
    
    if output_resize:
        total_width = int(total_width/output_resize_multiplier)
        max_height = int(max_height/output_resize_multiplier)
    
    combined_image.resize((total_width, max_height))
    
    if show:
        display(combined_image)
    # Save the combined image
    if outdir == None:
        outdir = indir

    if save:
        combined_image.save(outdir + marker + ".combined.png")
    
    # Close the input images
    image1.close()
    image2.close()
    image3.close()

### end of plots ###
### auxiliaries ###

def format_float(f, dec = 4):
    return f'{f:.{dec}f}'
def quantile_normalisation(df, col):
    return preprocessing.quantile_transform(df[col].values.reshape(-1, 1), output_distribution = 'normal', n_quantiles = df.shape[0])
def get_smallest_p(file):
    df = pd.read_csv(file, sep = ' ').dropna()
    df['p'] = df['p'].astype(float)
    df = df.sort_values(by = 'p').head(1)
    print(df)
    return df[['chr', 'pos']]
def encode_sex(s):
    if s == 'Female':
        return 2
    elif s == 'Male':
        return 1
    else:
        return 0
### end of auxiliaries ###

### phe ###

# Generate phe.txt file
def get_phenotype(df, col, quantile_transformation = True, convert_decimal = False, trailing_str = 'BOSON', save = False, save_path = None):
    df['Patient_ID'] = trailing_str + df['Patient_ID']
    df['IID'] = df['Patient_ID']
    if quantile_transformation:
        df[col] = quantile_normalisation(df, col)
    if convert_decimal:
        df[col] = df[col]/100
    phe = df.loc[:, ['Patient_ID', 'IID', col]]
    phe[col] = phe[col].apply(format_float)
    if save:
        phe.to_csv(save_path, header = None, index = None, sep = ' ')
    return None

# Plot original distribution and filter out outliers
def get_phe(col, df, save_name, plot = False, save_plot = False, save_path_plot = None, save_phe = False, save_path_phe = None):
    if save_plot and save_path_plot is None:
        save_path_plot = "graphs/facs-host/baseline/strip." + save_name + ".png"
    if save_phe and save_path_phe is None:
        save_path_phe = "boson_vcf/phenotypes/" + save_name + ".txt"
    tmp = sanity_check_filter4(df, col, plot = plot, save = save_plot, save_path = save_path_plot)
    get_phenotype(tmp, col, save = save_phe, save_path = save_path_phe)### end of phe ###
### parse plink GWASs results ###

def parse_plink(file, save = False, save_path = None):
    gwas = []
    for line in open(file): 
        lst = line.split('\n')[0].strip().split(' ')
        lstnew = [i for i in lst if i != '']
        gwas.append(lstnew)
    gwas.pop(0)
    gwas = pd.DataFrame(gwas, columns=['chr', 'ID', 'pos', 'alt', 'add', 'sample_size', 'beta', 'SE', 'p']).dropna()
    gwas = gwas[gwas['p'] != 'NA']
    gwas['p'] = gwas['p'].astype(float)
    gwas = gwas.sort_values(by = 'p', ascending = True).dropna()
    if save:
        gwas.to_csv(save_path, index = None, sep = ' ')
    return gwas
