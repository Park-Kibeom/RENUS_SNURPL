import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



class animate():
    def __init__(self, datagen, interval=10):
        """
        """
        self.time_interval = interval
        self.datagen = datagen
        self.fig, ax = plt.subplots()
        ax.set_title('Xenon oscillations', fontsize=20)

        self.line, = ax.plot([], [],"-r", linewidth = 2)  
        ax.set_ylim(0., 2.)
        ax.set_xlim(0, 1)
        self.time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        ax.grid(b=True)
        self.ax = ax

    def show(self):
        ani = animation.FuncAnimation(self.fig, self.refresh, self.datagen, interval=self.time_interval, blit=True, repeat=False, init_func=self.init)
        plt.show()

#save(filename, writer=None, fps=None, dpi=None, codec=None, bitrate=None, extra_args=None, metadata=None, extra_anim=None, savefig_kwargs=None)

    def init(self):
        self.line.set_data([], [])
        self.time_text.set_text('')
        return self.line, self.time_text

    def refresh(self,data):
        time, xdata, ydata = data
        self.line.set_xdata(xdata)
        self.line.set_ydata(ydata)  # update the data
        self.time_text.set_text('time = %.2fh'%(time))
        return self.line, self.time_text




#xdata = np.linspace(0,1,10)
#time = np.linspace(0,10,100)
#fun = lambda t,x: np.sin(np.pi*x) + 0.5*t*np.sin(2*np.pi*x) 
#ydata = np.array([fun(t,xdata) for t in time])
#animate(time,xdata,ydata)

