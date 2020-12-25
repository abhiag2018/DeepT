import sklearn
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_roc(saveto, name, labels, predictions):
    mpl.rcParams['figure.figsize'] = (12, 10)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fp, tp, _ = sklearn.metrics.roc_curve(labels, predictions)

    auc = sklearn.metrics.auc(fp, tp)
    plt.plot(100*fp, 100*tp, label = name+f"; auc : {auc:.2f}", linewidth=2, color=colors[0], linestyle='--')
    plt.legend()
    plt.xlabel('False positives [%]')
    plt.ylabel('True positives [%]')
    plt.xlim([0,100])
    plt.ylim([0,100])
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(saveto)
    plt.close('all')

def plot_prc(saveto, name, labels, predictions):
    mpl.rcParams['figure.figsize'] = (12, 10)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    pr, rec, _ = sklearn.metrics.precision_recall_curve(labels, predictions)

    auc = sklearn.metrics.auc(rec, pr)
    plt.plot(100*rec, 100*pr, label = name+f"; auc : {auc:.2f}", linewidth=2, color=colors[0], linestyle='--')
    plt.legend()
    plt.xlabel('Recall [%]')
    plt.ylabel('Precision [%]')
    plt.xlim([0,100])
    plt.ylim([0,100])
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(saveto)
    plt.close('all')