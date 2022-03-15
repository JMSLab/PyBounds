from matplotlib import pyplot as plt
import PyBounds.pyBounds as pyBounds


class Settings(pyBounds.Bounds):

    def plot_helper(self, x, lower_bounds, upper_bounds):
        fig, (ax, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={
            "width_ratios": [self.underline_B(), self.maxB - self.underline_B()]})
        ax.plot(x, lower_bounds, linewidth=.3)
        ax2.plot(x, lower_bounds, linewidth=.3)
        ax.plot(x, upper_bounds, linewidth=.3)
        ax2.plot(x, upper_bounds, linewidth=.3)

        ax.set_xlim(0, self.underline_B(), auto=True)
        ax2.set_xlim(self.underline_B(), self.maxB, auto=True)

        ax.spines['bottom'].set_linestyle((0, (8, 5)))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.yaxis.tick_left()
        ax2.yaxis.set_ticks_position('none')
        plt.subplots_adjust(wspace=0)

        fig.text(0.5, 0.04, 'Bound B', ha='center', va='center')
        ax.set_ylabel('Slope $\\theta$')
        plt.savefig('Bounds Interval Plot k=' + str(self.k) + ".png")
        plt.show()
