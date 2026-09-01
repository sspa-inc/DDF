#%%
import flopy
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from matplotlib import patheffects
from adjustText import adjust_text

#%% create MF grid
workspace = "mf"
print("Reading in spatial datasets")
wel = gpd.read_file("well.shp").sort_values("WellName")
riv = pd.read_csv("riv.csv")
riv = gpd.GeoDataFrame(riv, geometry=[LineString([(r.x1,r.y1), (r.x2, r.y2)]) for i,r in riv.iterrows()]).set_crs(epsg=26852)
grid = gpd.read_file("mf_grid.shp").set_crs(epsg=26852, allow_override=True)
gridpnt = grid.centroid
gridx = gridpnt.x
gridy = gridpnt.y
wel = gpd.read_file("well.shp").to_crs(epsg=26852).sort_values("WellName")
seg = gpd.read_file("seg.shp").to_crs(epsg=26852)
extent_img = [2435684.0683,2601365.7137,352035.2604,535699.0511]
extents = [
    (2471989,432771,2496222,452985),
    (2473029,392804,2518907,434347),
    (2471623,471373,2505355,498407),
    (2478930,446288,2508399,473322),
    (2493665,428022,2525326,455665),
    (2507528,444281,2557740,495396),
]
#%%
print('Contouring MODFLOW results and plotting drawdown')
levels=[0.5, 1, 2, 4, 9]
scalebar_left = [0.1, 0.1, 0.6, 0.7, 0.6, 0.7]
scalebar_size = np.array([2, 4, 3, 2, 3, 4])

fig, axs = plt.subplots(2, 3, figsize=(6.75, 4.70), tight_layout=True)
for i, ax in enumerate(axs.flatten()):
    gg = gpd.sjoin(grid, gpd.GeoDataFrame(geometry=wel.iloc[i:i+1].buffer(33123)))
    gg.plot(ax=ax, facecolor="none", edgecolor="grey", linewidth=0.25, )
    seg.plot(ax=ax, color="b")
    wel.iloc[i:i+1].plot(ax=ax, color="r", markersize=8)
    contour = gpd.read_file(f"well{i+1}_cc_18250.000.geojson", )
    contour=contour.set_crs(epsg=26852,allow_override=True)
    contour[contour.Level.isin(levels)].plot(ax=ax, color="r", linewidth=1.5, zorder=10, alpha=0.5)
    head = -flopy.utils.HeadFile(f"mf_wel{i+1}/Output.hds", ).get_data(idx=0)[0,0]
    ax.tricontour(gridx, gridy, head, levels=levels, linewidths=1.0, colors='k', linestyles="--", zorder=20)
    extent = extents[i]
    dx = extent[2] - extent[0]
    dy = extent[3] - extent[1]
    dx = max(dx, dy)
    dy = dx
    ax.text(0.05 if i!=1 else 0.75, 0.95, f"Well {i+1}", va="top", transform=ax.transAxes, fontsize="medium",
            path_effects=[patheffects.withStroke(linewidth=3, foreground="w")])


    ax.text(extent[0]+dx*scalebar_left[i]+3280.8*scalebar_size[i]/2, extent[1]+dx*0.055, f"{scalebar_size[i]} km", fontsize="small", # + ("s" if scalebar_size[i]>1 else ""),
            ha="center", va="bottom", path_effects=[patheffects.withStroke(linewidth=3, foreground="w")])
    ax.plot([extent[0]+dx*scalebar_left[i], extent[0]+dx*scalebar_left[i]+3280.8*scalebar_size[i]], [extent[1]+dx*0.05,]*2, color="k",
            path_effects=[patheffects.withStroke(linewidth=3, foreground="w")])
    ax.set(xlim=(extent[0], extent[0]+dx), ylim=(extent[1], extent[1]+dx))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

ax = axs[0,0]
l1 = ax.plot([0],[0], color="r", linewidth=2, zorder=10, alpha=0.5, label="DDF")[0]
l2 = ax.plot([0],[0], linewidth=1.0, color='k', linestyle="--", label="MODFLOW")[0]
ax.legend(handles=[l1, l2], fontsize="small", )
fig.tight_layout(pad=0.2)
fig.savefig("plot_ddf_dd.png", dpi=600, )
#%%
print("Plotting MODFLOW grid")
extent = [2465424.41107268,2547735.15914244,397071.49541570,508563.42936919]
fig, ax = plt.subplots(figsize=(3.248, 4), tight_layout=True)
seg.plot(ax=ax, color="b", zorder=2)
wel.plot(ax=ax, color="r", markersize=10, zorder=3)
labels = []
for i, r in wel.iterrows():
    labels.append(ax.text(r.geometry.x, r.geometry.y, r.WellName, zorder=4, fontsize="x-small", path_effects=[patheffects.withStroke(linewidth=2, foreground="w")]))
adjust_text(labels)
grid.plot(ax=ax, facecolor="none", edgecolor="grey", linewidth=0.25, zorder=1)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

ax.text(2.538e6,498000, u'\u25B2\nN', fontsize="medium", zorder=4, path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])
ax.text(2477074+3280.8*5, 400100+500, "10 km", zorder=4, fontsize="medium", # + ("s" if scalebar_size[i]>1 else ""),
        ha="center", va="bottom", path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])
ax.plot([2477074, 2477074+3280.8*10], [400100,]*2, color="k", zorder=4,
        path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])
ax.set(xlim=(extent[0], extent[1]), ylim=(extent[2], extent[3]))
fig.savefig("plot_map_grid.png", dpi=600, )


#%%
print('Plotting stream depletion factor at Year 50')
plt.style.use("seaborn-v0_8-whitegrid")
sds = {}
for i in range(6):
    sd = []
    ddf = pd.read_csv(f"well{i+1}_dl.csv", )
    sd.append(ddf.loc[:, '18250.0000'].sum()/20000)

    lst = flopy.utils.Mf6ListBudget(f"mf_wel{i+1}/Output.lst")
    df1, df2 = lst.get_dataframes(diff=True)
    sd.append(df2.chd.iloc[-1]/np.abs(df2.wel.iloc[-1]))

    sds[f"Well{i+1}"] = sd

fig, ax = plt.subplots(figsize=(8, 4))
ax = (pd.DataFrame(sds, index=["DDF", "MODFLOW"])*100).T.plot.bar(ax=ax)
ax.set(ylabel="Stream depletion factor (%)", )
ax.figure.savefig("plot_dl.png", dpi=600, bbox_inches="tight")
pd.DataFrame(sds, index=["DDF", "MODFLOW"]).to_csv("plot_dl.csv")
#%%
print('Plotting stream depletion timeseries')
fig, ax = plt.subplots(figsize=(8, 5))
plt.style.use("seaborn-v0_8-whitegrid")

for i in range(6):
    ax.clear()
    ddf = pd.read_csv(f"well{i+1}_dl.csv", ).iloc[:,3:].sum()
    ddf.index = pd.to_numeric(ddf.index) / 365
    ddf = ddf / 20000 * 100
    ddf[(ddf>0)&(ddf<1e6)].plot(ax=ax, label="DDF")

    lst = flopy.utils.Mf6ListBudget(f"mf_wel{i+1}/Output.lst")
    df1, df2 = lst.get_dataframes(diff=True)
    ax.plot(np.array(lst.get_times())/365, df2.chd/df2.wel.abs()*100, label="MODFLOW")

    ax.legend()
    ax.set(xlabel="Elasped time (years)", ylabel="Stream depletion factor (%)")
    fig.savefig(f"well{i+1}_dl.png", dpi=600, bbox_inches="tight")
