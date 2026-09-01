import geopandas as gpd
import pandas as pd
import numpy as np

#%%
wel = gpd.read_file("well.shp").sort_values("WellName")
xy_wel = np.array([wel.geometry.x, wel.geometry.y]).T.astype(int)

wel_file = pd.DataFrame(xy_wel, columns=["x", "y"])
wel_file["Transmissivity"] = 1000
wel_file["Storativity"] = 0.1
wel_file["Rate"] = 20000
wel_file["Timeoff"] = 18250

for i in range(6):
    print(f"Creating files for well {i+1}")
    w = wel_file.loc[i:i]
    w.to_csv(f"well{i+1}.csv", index=False)

    with open(f"well{i+1}.config", "w") as f:
        f.write("1.0 1e-3 1e-4 0.5 0\n")         # lfac, ddriv_min
        f.write("50 1.1 0\n")         # cellsize1, cellmutl, gridrmax
        f.write("0.5\n")              # contourlevel
        #f.write("-1 \n")              # the contour level; not used if grid file is provided
        f.write(f"well{i+1}.csv\n")   # rite well info into this file; pumping is positive for rates
        f.write("riv.csv\n")          # river reach file
        f.write("times.csv\n")        # write the time when the drawdown contour will be created
        f.write(f"well{i+1}\n")       # prefix for output files