import matplotlib.pyplot as plt

plt.figure(num=None, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')
legends = []

def lstplot(xdata,ydata,xlabel="x",ylabel="y",legend="?",style="r-"):
    """
    """
    plt.plot(xdata,ydata,style, linewidth = 2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(b=True)
    legends.append(legend)
    plt.legend(legends, loc='best')
