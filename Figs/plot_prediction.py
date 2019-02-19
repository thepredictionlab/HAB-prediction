import numpy as np
import matplotlib.pylab as plt
from matplotlib.pyplot import figure, show, rc
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from collections import OrderedDict
import matplotlib.colors as mcolors
import seaborn as sns


####### PREDICTION PLOT
sns.set_style('white')
sns.set_context('notebook')
sns.set_style("whitegrid", {'axes.grid' : False})

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])

    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

c = mcolors.ColorConverter().to_rgb
rvb = make_colormap(
[c('red'), 0.125, c('red'), c('orange'), 0.25, c('orange'),c('green'),0.5, c('green'),0.7, c('green'), c('blue'), 0.75, c('blue')])

#! Plot
fig, ax = plt.subplots(figsize=(8,5))
Cert = np.asarray([.1,.6,.85,.4,.1])
Prob = np.asarray([.6,.2,.1,.05,.05])
objects = ('Worse alot', 'Worse a little', 'Same', 'Better a little', 'Better alot')
N = len(objects)
y_pos = len(objects) - np.arange(len(objects)) - 1
colors = rvb(1-Cert)

# Add colormap
y = np.unique(Cert)
newcmp = ListedColormap(rvb(1-y))
plot = plt.scatter(y, y, c = y, cmap = newcmp,alpha=0.5)
plt.clf()
plt.colorbar(plot)

#! Barplot
plt.barh(y_pos, Prob, align='center', alpha=0.5,color=colors)
plt.yticks(y_pos, objects)
plt.xlabel('Probability')
plt.title('HAB Condition Forecast')

# Turns off grid on the left Axis.
ax.grid(False)

#! Add current date
sns.despine(ax=ax, offset=10)
fig.tight_layout()
plt.show()
plt.savefig("./PNG/Fig_prediction.png",dpi=600)

