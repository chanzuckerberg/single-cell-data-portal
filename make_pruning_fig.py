import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import numpy as np

rng = np.random.default_rng(7)

READ = "#2E7D32"; READ_FILL = "#D9EAD3"
SKIP = "#B5B5B5"; SKIP_FILL = "#F1F1F1"
BAND = "#4A78B5"; BLOOM = "#C77D0A"; BLOOM_FILL = "#FDE9C8"

fig = plt.figure(figsize=(15, 8.4))
gs = fig.add_gridspec(2, 2, height_ratios=[6.2, 1.25], hspace=0.30, wspace=0.16)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axT = fig.add_subplot(gs[1, :]); axT.axis("off")

fig.suptitle("Why the same selective lookup prunes differently: TileDB vs ClickHouse",
             fontsize=15, fontweight="bold", y=0.99)

# ---------------- Panel A: TileDB ----------------
axA.set_title("TileDB — coordinate-space pruning\nR-tree over tile bounding boxes (~10k cells/tile)",
              fontsize=12.5, fontweight="bold")
axA.set_xlim(0, 10); axA.set_ylim(0, 10)
axA.set_xlabel("gene  →", fontsize=11); axA.set_ylabel("tissue  →", fontsize=11)
axA.set_xticks([]); axA.set_yticks([])

gq = (3.1, 5.3); tq = (4.4, 7.1)   # query: gene IN(...) AND tissue IN(...)
axA.add_patch(Rectangle((gq[0], 0), gq[1]-gq[0], 10, color=BAND, alpha=0.09, zorder=0))
axA.add_patch(Rectangle((0, tq[0]), 10, tq[1]-tq[0], color=BAND, alpha=0.09, zorder=0))
axA.add_patch(Rectangle((gq[0], tq[0]), gq[1]-gq[0], tq[1]-tq[0], color=BAND, alpha=0.20, zorder=1))
axA.text(gq[0]+ (gq[1]-gq[0])/2, 0.35, "gene IN(…)", ha="center", fontsize=9, color=BAND)
axA.text(0.25, tq[0]+(tq[1]-tq[0])/2, "tissue\nIN(…)", va="center", fontsize=9, color=BAND)

tiles = [(1.6, 2.0, 1.7, 1.7), (4.3, 5.9, 1.9, 1.9), (2.1, 7.7, 1.6, 1.5),
         (6.7, 2.4, 1.8, 1.7), (7.7, 7.2, 1.7, 1.7), (4.7, 3.0, 1.6, 1.6)]
def overlaps(t):
    cx, cy, w, h = t; x0, x1 = cx-w/2, cx+w/2; y0, y1 = cy-h/2, cy+h/2
    return not (x1 < gq[0] or x0 > gq[1] or y1 < tq[0] or y0 > tq[1])
for t in tiles:
    cx, cy, w, h = t; x0, y0 = cx-w/2, cy-h/2; hit = overlaps(t)
    pts = rng.uniform([x0+0.18, y0+0.18], [x0+w-0.18, y0+h-0.18], size=(16, 2))
    axA.scatter(pts[:, 0], pts[:, 1], s=6, color=(READ if hit else SKIP), alpha=0.85, zorder=3)
    axA.add_patch(Rectangle((x0, y0), w, h, fill=False, lw=(2.4 if hit else 1.2),
                  edgecolor=(READ if hit else SKIP), linestyle=("-" if hit else (0, (4, 3))), zorder=4))
axA.text(4.3, 5.9, "READ", ha="center", va="center", fontsize=8.5, fontweight="bold", color=READ, zorder=5)
axA.text(4.7, 3.0, "READ", ha="center", va="center", fontsize=8.5, fontweight="bold", color=READ, zorder=5)
axA.text(9.6, 9.4, "MBR bounds\nALL dims per tile", ha="right", va="top", fontsize=8.5, color="#555")
axA.text(5.0, -0.9, "Intersects predicates on several dims at once → dimension order doesn't matter.",
         ha="center", fontsize=9.5, color="#333")

# ---------------- Panel B: ClickHouse ----------------
axB.set_title("ClickHouse — sorted-prefix pruning\nsparse primary index over rows sorted by the key",
              fontsize=12.5, fontweight="bold")
axB.set_xlim(0, 15); axB.set_ylim(1.4, 10.4); axB.axis("off")

n = 14
read_idx = {5, 6, 7}
# primary-key strip (on-prefix gene lookup)
axB.text(0.3, 9.9, "ORDER BY (organism, gene, tissue)  —  rows in 8192-row granules",
         fontsize=9.5, color="#333", fontweight="bold")
y0 = 6.2
# sparse index marks (above strip)
for i in range(n+1):
    axB.plot([0.3+i, 0.3+i], [y0+1.5, y0+1.85], color="#888", lw=0.8)
axB.text(0.3, y0+2.0, "sparse index marks (granule boundary key values, kept in memory)",
         fontsize=8.2, color="#666")
for i in range(n):
    read = i in read_idx
    axB.add_patch(Rectangle((0.3+i, y0), 0.92, 1.5,
                  facecolor=(READ_FILL if read else SKIP_FILL),
                  edgecolor=(READ if read else SKIP), lw=(2.2 if read else 1.0)))
axB.annotate("gene IN(…): binary-search marks → read 3 / 42,969 granules",
             xy=(0.3+6.5, y0), xytext=(0.3+6.5, y0-1.0), ha="center", fontsize=9.5, color=READ,
             arrowprops=dict(arrowstyle="-", color=READ, lw=1.2))

# off-prefix strip (cell_type) — bloom skip index, coarse/probabilistic
y1 = 2.9
bloom_idx = {1, 4, 5, 8, 11}
for i in range(n):
    hit = i in bloom_idx
    axB.add_patch(Rectangle((0.3+i, y1), 0.92, 1.5,
                  facecolor=(BLOOM_FILL if hit else SKIP_FILL),
                  edgecolor=(BLOOM if hit else SKIP), lw=(1.8 if hit else 1.0),
                  linestyle=("-" if hit else "-")))
axB.text(0.2, y1+1.75, "off-prefix filter: cell_type  —  primary index can't help",
         fontsize=9.5, color="#333")
axB.text(0.3+7, y1-0.85,
         "falls back to a data-skipping (bloom) index → coarse, probabilistic  (~1.2× slower)",
         ha="center", fontsize=9.5, color=BLOOM)

# ---------------- Takeaway ----------------
axT.set_xlim(0, 1); axT.set_ylim(0, 1)
axT.add_patch(Rectangle((0.005, 0.05), 0.99, 0.9, facecolor="#F7F7F7", edgecolor="#DDD"))
axT.text(0.02, 0.72,
         "Takeaway:  TileDB gives flexible multi-dimensional slicing natively (any dimension order).  "
         "ClickHouse optimizes one favored access path — the sort-key prefix — with coarser skip-indexes for everything else.",
         fontsize=11, color="#111", fontweight="bold", va="center")
axT.text(0.02, 0.30,
         "A well-chosen sort key mimics coordinate tiling by clustering data contiguously; a poor one triggers scans.  "
         "Lance failed a third way: its scalar index returned row-ids, then did a scattered take over bitmap space with no contiguous "
         "clustering — a ~35 ms floor.  The deciding property is the on-disk LAYOUT, not the “columnar” label.",
         fontsize=10, color="#333", va="center")

fig.savefig("/tmp/pruning_fig.png", dpi=150, bbox_inches="tight", facecolor="white")
fig.savefig("/tmp/pruning_fig.svg", bbox_inches="tight", facecolor="white")
print("ok")
