import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

padding = 1.5
leftMargin = 0.1
rightMargin = 0.1
upMargin = 0.5
bottomMargin = 0.5
def adjustLegendPadding(plt,ax, scaling = 1.5):
    plt.tight_layout(pad = padding * scaling)    
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.7, chartBox.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1, prop = {"size":6})        
    params = {'legend.fontsize': 5,
            'legend.handlelength': 2}
    plt.rcParams.update(params)

def angularPowerAxisLabel(plt,ax, scaled = False):
    """
        used to scale the axis

    Parameters:
    -----------------
    
    scaled: boolean
        default False, if True, will extra times 4 Pi / N^2
    """
    plt.xlabel("$\ell$")
    if (scaled == False):
        plt.ylabel("$(2\ell+1)C_{\ell}$")
    else:
        # plt.ylabel(r'$\frac{4\pi(2\ell+1)C_{\ell}}{N^2}')
        # plt.tight_layout(pad = padding * 2)
        plt.ylabel(r"$\frac{4 \pi(2\ell+1)C_{\ell}}{N^2}$", labelpad = -1)

    plt.grid()    

def correlationLabel(plt,ax):
    print("currently not supported")
    pass
