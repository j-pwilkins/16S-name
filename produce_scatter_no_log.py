import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main():
    data_df = pd.read_csv('Cs_Decimal_Matrix.csv')
    produce_scatter_plot(data_df)

def produce_scatter_plot(data_df):
    plt.figure(figsize=(12, 10))
    plot_region(data_df, 'Full', 'V1V3', 'blue', 1, 'V1V3 Mapped Scores vs Full 16S Mapped Scores',
                'Full 16S Mapped Scores', 'V1V3 Mapped Scores', (0, 0.31), (0, 0.31), plot_color="b-")
    plot_region(data_df, 'Full', 'V3V4', 'green', 2, 'V3V4 Mapped Scores vs Full 16S Mapped Scores',
                'Full 16S Mapped Scores', 'V3V4 Mapped Scores', (0, 0.31), (0, 0.31), plot_color="g-")
    plot_region(data_df, 'Full', 'V4', 'red', 3, 'V4 Mapped Scores vs Full 16S Mapped Scores',
                'Full 16S Mapped Scores', 'V4 Mapped Scores', (0, 0.31), (0, 0.31), plot_color="r--")
    plot_trendlines(data_df)
    plt.tight_layout()
    plt.show()

def plot_region(df, x, y, color, subplot_index, title, xlabel, ylabel, xlim, ylim, plot_line=True, linewidth=0.8, dashes=[2, 5], plot_color=None):
    plt.subplot(2, 2, subplot_index)
    plt.scatter(df[x], df[y], color=color, s=2)
    plt.title(title)
    plt.xlabel(xlabel, color='orange')
    plt.ylabel(ylabel, color=color, fontweight='bold')
    plt.xlim(xlim)
    plt.ylim(ylim)
    if plot_line:
        z = np.polyfit(df[x], df[y], 1)
        p = np.poly1d(z)
        if plot_color is not None:
            plt.plot(df[x], p(df[x]), plot_color, linewidth=linewidth, dashes=dashes)
        else:
            plt.plot(df[x], p(df[x]), linewidth=linewidth, dashes=dashes)

def plot_trendlines(data_df):
    plt.subplot(2, 2, 4)
    z = np.polyfit(data_df['Full'], data_df['V1V3'], 1)
    p = np.poly1d(z)
    plt.plot(data_df['Full'], p(data_df['Full']), "b--", linewidth=0.8, dashes=[1, 5])
    z = np.polyfit(data_df['Full'], data_df['V3V4'], 1)
    p = np.poly1d(z)
    plt.plot(data_df['Full'], p(data_df['Full']), 'g--', linewidth=0.8, dashes=[1, 5])
    z = np.polyfit(data_df['Full'], data_df['V4'], 1)
    p = np.poly1d(z)
    plt.plot(data_df['Full'], p(data_df['Full']), 'r--', linewidth=0.8, dashes=[1, 5])
    plt.xlim(0, 0.31)
    plt.ylim(0, 0.31)
    plt.xlabel('Full 16S Mapped Scores', color='orange')
    plt.ylabel('Regional Mapped Scores', color='black', fontweight='bold')
    plt.title('Regional Trendlines for Mapped Scores')

## helper functions
def read_csv(filename):
    try:
        df = pd.read_csv(filename)
        return df
    except FileNotFoundError:
        print(
            "The file '{}' could not be found. Make sure the file is in the correct location and try again.".format(
                filename))
        exit()
    except pd.errors.EmptyDataError:
        print("The file '{}' is empty. Make sure the file contains data and try again.".format(filename))
        exit()
    except:
        print("An unknown error occurred while trying to read the file '{}'.".format(filename))
        exit()

def write_csv(df, filename):
    df.to_csv(filename, sep=',', index=False)

if __name__ == '__main__':
    main()