import sklearn
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_roc(saveto, name, labels, predictions, **kwargs):
    mpl.rcParams['figure.figsize'] = (12, 10)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fp, tp, _ = sklearn.metrics.roc_curve(labels, predictions)

    plt.plot(100*fp, 100*tp, label=name, linewidth=2, color=colors[0], linestyle='--', **kwargs)
    plt.xlabel('False positives [%]')
    plt.ylabel('True positives [%]')
    plt.xlim([0,100])
    plt.ylim([0,100])
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(saveto)
    plt.close('all')