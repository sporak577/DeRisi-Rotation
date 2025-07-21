import pandas as pd
import numpy as np
import pycountry
import plotly.express as px

# 1. Load & clean
df = pd.read_csv('file.csv') #update with your toy/real dataframe here
df["geo_loc_name"] = df["geo_loc_name"].replace("nan", np.nan)

# 2. Ask user which column to facet by, and which values to include
print("Columns you can facet on:", ", ".join(df.columns))
facet_col = input("Enter the column to facet by: ").strip()
if facet_col not in df.columns:
    raise KeyError(f"Column '{facet_col}' not found in your data.")

raw_vals = input(f"Enter the values from '{facet_col}' you want (comma-separated): ")
facet_vals = [v.strip() for v in raw_vals.split(",") if v.strip()]

# 3. Subset to just those values
df = df[df[facet_col].isin(facet_vals)]
print(f"\nKept {len(df)} rows matching {facet_col} ∈ {facet_vals}\n")

# 4. ISO-3 conversion
def country_to_iso3(name):
    if pd.isna(name):
        return None
    # Strip anything after colon, comma, or dash
    name_clean = str(name).split(":")[0].split(",")[0].split("-")[0].strip()
    try:
        return pycountry.countries.lookup(name_clean).alpha_3
    except LookupError:
        return None
    
df["iso_alpha"] = df["geo_loc_name"].apply(country_to_iso3)
df = df.dropna(subset=["iso_alpha"])

# 5. Aggregate counts by country & facet value
agg = (
    df
    .groupby(["iso_alpha", facet_col], as_index=False)["count"]
    .sum()
)

# 6. Draw faceted choropleth
fig = px.choropleth(
    agg,
    locations="iso_alpha",
    color="count",
    facet_col=facet_col,
    projection="natural earth",
    title=f"Count by Country, faceted by {facet_col}",
    hover_name="iso_alpha",
    hover_data={"count": True},
    color_continuous_scale="Viridis",
)
fig.update_geos(showframe=False, showcoastlines=True)
fig.show()
