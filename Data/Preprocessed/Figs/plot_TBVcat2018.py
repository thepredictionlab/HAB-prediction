import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# replace with whatever path to the 2018results file
dat = np.load('./Results/TBVcat/2018Pred/2018results.npz')
post = dat['POSTERIOR']
bins = dat['BINS']
date = dat['DATE']
obs = dat['OBS']

m = []
u = []
l = []

for k in np.arange(279,290):
	data = post.flatten()[0][k]
	cdf = np.cumsum(data)
	a1 = np.max(np.where(cdf<=0.5))
	a2 = np.min(np.where(cdf>0.5))
	y1 = cdf[a1]
	y2 = cdf[a2]
	x1 = bins[a1]
	x2 = bins[a2]
	sl = (x2-x1)/(y2-y1)
	m.append(sl*(0.5-y1) + x1)
	a1 = np.max(np.where(cdf<0.05))
	a2 = np.min(np.where(cdf>=0.05))
	y1 = cdf[a1]
	y2 = cdf[a2]
	x1 = bins[a1]
	x2 = bins[a2]
	sl = (x2-x1)/(y2-y1)
	l.append(sl*(0.05-y1) + x1)
	a1 = np.max(np.where(cdf<=0.95))
	a2 = np.min(np.where(cdf>0.95))
	y1 = cdf[a1]
	y2 = cdf[a2]
	x1 = bins[a1]
	x2 = bins[a2]
	sl = (x2-x1)/(y2-y1)
	u.append(sl*(0.95-y1) + x1)

plt.plot(date[279:290],m)
plt.plot(date[279:290],u)
plt.plot(date[279:290],l)
plt.scatter(date[279:290],obs[279:290])
plt.ylim(-1.5,1.5)
plt.title('Predicted and observed percentage change \n in cyanobacteria biovolume during 2018 season.')
plt.show()


