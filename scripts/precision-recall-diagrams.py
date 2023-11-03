import numpy as np
import pylab as pl
import itertools as it

def fmeasure(p, r):
    """Calculates the fmeasure for precision p and recall r."""
    return 2*p*r / (p+r)

def _fmeasureCurve(f, p):
    """For a given f1 value and precision get the recall value."""
    return f * p / (2 * p - f) if (2 * p - f) != 0 else 0

def _plotFMeasures(fstepsize=.1, stepsize=0.001):
    """Plots 10 fmeasure Curves into the current canvas."""
    p = np.arange(0., 1., stepsize)[1:]
    for f in np.arange(0., 1., fstepsize)[1:]:
        points = [(x, _fmeasureCurve(f, x)) for x in p if 0 < _fmeasureCurve(f, x) <= 1.5]
        xs, ys = zip(*points)
        curve, = pl.plot(xs, ys, "--", color="gray", linewidth=0.5)
        pl.annotate(r"$f=%.1f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.05, ys[-10] - 0.035), size="small", color="gray")

colors = "bgrcmyk"
markers = "so^>v<dph8"

def plotPrecisionRecallDiagram(title="title", points=None, labels=None, loc="center right"):
    if labels:
        ax = pl.axes([0.1, 0.1, 0.7, 0.8])
    else:
        ax = pl.gca()
    pl.title(title)
    pl.xlabel("Precision")
    pl.ylabel("Recall")
    _plotFMeasures()

    if points is not None and len(points) > 0:
        getColor = it.cycle(colors).__next__
        getMarker = it.cycle(markers).__next__

        scps = []
        for i, (x, y) in enumerate(points):
            label = None
            if labels:
                label = labels[i]
            scp = ax.scatter(x, y, label=label, s=50, linewidths=0.75, facecolor=getColor(), alpha=0.75, marker=getMarker())
            scps.append(scp)

        if labels:
            pl.legend(scps, labels, loc=(1.01, 0), scatterpoints=1, numpoints=1, fancybox=True)
    pl.axis([-0.02, 1.02, -0.02, 1.02])


prs = np.array([
    [0.7151839776432231, 0.3629165681871898],  # Bulk ONT
    [0.5886442641946698, 0.8704797920113448],  # Bulk Pacbio
    [0.25528346192552837, 0.22482864571023398]  # MDA-ONT
])

labels = ["Bulk ONT", "Bulk Pacbio", "MDA-ONT"]

pl.figure(figsize=(10, 6))

plotPrecisionRecallDiagram("Benchmarking against PC3 ground truth", prs, labels)
pl.savefig('result_plot.pdf')
pl.show()


