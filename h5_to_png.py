#!/usr/bin/env python
import numpy as np
import mapio
import argparse

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import scipy.ndimage.morphology as morph

def prep_matrix(a):
    np.fill_diagonal(a, 1.)
    d = morph.grey_closing(a, size=(3,3))
    m = morph.binary_erosion(a, structure=np.ones((3,3)))
    m = morph.binary_propagation(m, mask=a).astype(np.float)
    np.fill_diagonal(m, 0)
    a *= m


parser = argparse.ArgumentParser('H5 contact map to PNG image')
parser.add_argument('--dual', action='store_true', default=False, help='Plot dual raw/filter side by side')
parser.add_argument('--filter', action='store_true', default=False, help='Apply enchancement filter')
parser.add_argument('--dpi', type=int, default=300, help='Resolution in DPI')
parser.add_argument('--edge-locs', help='Location of cluster borders')
parser.add_argument('h5_in', help='Input H5 contact map')
parser.add_argument('out_name', help='Output file name')
args = parser.parse_args()


cm, names = mapio.read_map(args.h5_in, fmt='h5', has_names=True)
cm = cm.todense()
print "read {0}x{0} map".format(len(cm))
print 'weight ', np.sum(cm)

m = [cm]
if args.dual:
    cm2 = cm.copy()
    prep_matrix(cm2)
    m.append(cm2)
elif args.filter:
    prep_matrix(cm)


plt.close()
plt.style.use('ggplot')
fig, axes = plt.subplots(1, len(m))
if len(m) == 1:
    axes = [axes]

fig.set_figheight(14)
fig.set_figwidth(14)


if args.edge_locs:
    locs = []
    with open(args.edge_locs, 'r') as in_h:
        for line in in_h:
            line = line.strip()
            if not line:
                break
            locs.append(int(line))
    ticks = ticker.FixedLocator(locs)

majlabels = ticker.FixedFormatter('')

plt.yticks(fontsize=10)
plt.xticks(fontsize=12, rotation=45)
plt.tick_params(axis='both', which='minor', labelsize=14)

for i, axi in enumerate(axes):
    
    axi.yaxis.set_major_formatter(majlabels)

    if args.edge_locs:
        axi.xaxis.set_major_locator(ticks)
        axi.yaxis.set_major_locator(ticks)

    axi.grid(linewidth=.1, color='0.2', linestyle=(0, (5,5)))
    axi.imshow(np.log(m[i] + .1), interpolation='none', origin='lower', cmap=matplotlib.cm.gnuplot2_r)
    
plt.tight_layout()
plt.savefig(args.out_name, dpi=args.dpi)

