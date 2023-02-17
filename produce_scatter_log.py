import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main():
    data_df = pd.read_csv('Cs_Decimal_Matrix.csv')
    produce_scatter_plot(data_df)

def produce_scatter_plot(data_df):
    plt.figure(figsize=(12, 10))
    # These are the limits selected for the non-log version. Chosen as 0.3 max value and leaves room from top of chart
    xlim = (0, 0.31)
    ylim = (0, 0.31)
    plot_region(data_df, 'Full', 'V1V3', 'blue', 1, 'V1V3 Mapped Scores vs Full 16S Mapped Scores',
                'Full 16S Mapped Scores', 'V1V3 Mapped Scores', xlim, ylim, plot_color="b-")
    plot_region(data_df, 'Full', 'V3V4', 'green', 2, 'V3V4 Mapped Scores vs Full 16S Mapped Scores',
                'Full 16S Mapped Scores', 'V3V4 Mapped Scores', xlim, ylim, plot_color="g-")
    plot_region(data_df, 'Full', 'V4', 'red', 3, 'V4 Mapped Scores vs Full 16S Mapped Scores',
                'Full 16S Mapped Scores', 'V4 Mapped Scores', xlim, ylim, plot_color="r--")
    plot_trendlines(data_df, xlim, ylim)
    plt.tight_layout()
    plt.show()

def plot_region(df, x, y, color, subplot_index, title, xlabel, ylabel, xlim, ylim, plot_color=None):
    plt.subplot(2, 2, subplot_index)
    plt.scatter(df[x], df[y], color=color, s=2)
    plt.title(title)
    plt.xlabel(xlabel, color='orange')
    plt.ylabel(ylabel, color=color, fontweight='bold')
    # plots the log scales onto the axes
    plt.xscale("log")
    plt.yscale("log")
    # set limits of chart - I don't know how this works, but without these lines plot is confined to portion of plot area
    plt.xlim(xlim)
    plt.ylim(ylim)
    # transforms the data. 0.00000001 is added to each value because there are zeroes which otherwise break the code
    x_log = np.log10(df[x] + 0.00000001)
    y_log = np.log10(df[y] + 0.00000001)
    # fit the trendlines
    z = np.polyfit(x_log, y_log, 1)
    p = np.poly1d(z)
    x_fit = np.linspace(x_log.min(), x_log.max(), 50)
    y_fit = 10 ** p(x_fit)
    plt.plot(10 ** x_fit, y_fit, plot_color, linewidth=0.8, dashes=[2, 5])
    add_clade_lines(0.05, xlim, ylim, ls='dotted', linewidth=1, color='black', label='species')
    add_clade_lines(0.20, xlim, ylim, ls='dotted', linewidth=1, color='black', label='genus')

def add_clade_lines(pos, xlim, ylim, ls, linewidth, color, label):
    plt.plot([pos, pos], [ylim[0], pos], ls=ls, linewidth=linewidth, color=color)
    plt.plot([xlim[0], pos], [pos, pos], ls=ls, linewidth=linewidth, color=color)
    plt.text(pos, pos/10, label, ha='center', va='bottom', fontsize=8, color=color)

    # I don't need the points plotted as below, but leaving code as I might want to plot some slightly different things
    # plt.plot([pos], [pos], marker='o', markersize=4, color='black')


def plot_trendlines(df, xlim, ylim):
    plt.subplot(2, 2, 4)
    plt.scatter(df['Full'], df['V1V3'], color='blue', s=0)
    plt.scatter(df['Full'], df['V3V4'], color='green', s=0)
    plt.scatter(df['Full'], df['V4'], color='red', s=0)
    plt.title('Regional Trendlines for Mapped Scores')
    plt.xlabel('Full 16S Mapped Scores', color='orange')
    plt.ylabel('Regional Mapped Scores', color='black', fontweight='bold')
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plot_trendline(df['Full'], df['V1V3'], "b")
    plot_trendline(df['Full'], df['V3V4'], "g")
    plot_trendline(df['Full'], df['V4'], "r")


def plot_trendline(x, y, color):
    x_log = np.log10(x + 0.00000001)
    y_log = np.log10(y + 0.00000001)
    z = np.polyfit(x_log, y_log, 1)
    p = np.poly1d(z)
    x_fit = np.linspace(x_log.min(), x_log.max(), 50)
    y_fit = 10 ** p(x_fit)
    plt.plot(10 ** x_fit, y_fit, color + '--', linewidth=0.8, dashes=[2, 5])


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