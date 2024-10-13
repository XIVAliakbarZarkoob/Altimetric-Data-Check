import requests
import os
import json

# ---- Niger -----------
# REGION_LAT = (2, 20)
# REGION_LON = (-15, 15)

# --- Ganges-Brahmaputra -----
REGION_LAT = (20, 33)
REGION_LON = (72, 100)


file_path = 'manifest_clms_global_wl_rivers_v2_daily_geojson_latest.txt'
save_directory = 'downloads'
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

# with open(file_path, "r") as file:
#     links = file.readlines()
    
for link in links:
    url = link.strip()
    try:
        response = requests.get(url)
        response.raise_for_status()  
        coord = json.loads(response.content)['geometry']['coordinates']
        if REGION_LON[0]<coord[0]<REGION_LON[1] and REGION_LAT[0]<coord[1]<REGION_LAT[1]:
            filename = url.split("/")[-1]  
            with open(f"{save_directory}/{filename}", "wb") as file:
                file.write(response.content)
            print(f"Downloaded: {filename}")
    except Exception as e:
        print(f"Failed to download {url}: {e}")