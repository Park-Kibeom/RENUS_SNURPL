import matplotlib.pyplot as plt

class lstplot:
    def __init__(self):
        self.plt = plt
        self.plt.figure(num=None, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')

    def addline(self,xdata,ydata,xlabel="x",ylabel="y",legend="?",style="r-"):
        """
        """
        self.plt.plot(xdata,ydata,style, linewidth = 2)
        self.plt.xlabel(xlabel)
        self.plt.ylabel(ylabel)
        self.plt.grid(b=True)
    #    legends.append(legend)
    #    plt.legend(legends, loc='best')
