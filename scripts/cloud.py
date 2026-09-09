from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import gc
import geopandas as gpd
import numpy as np
import omnicloudmask
import pystac_client
import rasterio
import torch

from affine import Affine
from rasterio.enums import Resampling
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union
from shapely.validation import make_valid

torch.set_num_threads(16)

tiles = [
'030008',
'026009',
'025009',
'008015',
'025008',
'007015',
'005014',
'006017',
'003015',
'002014',
'006014',
'013011',
'004016',
'005016',
'006016',
'009013',
'003016',
'010013',
'008013'
]

datetime_range = "2026-08-13/2026-08-15"

def read_band(href: str, resolution: float = 20.0):

    with rasterio.open(href) as src:

        scale = src.res[0] / resolution

        out_shape = (
            int(src.height * scale),
            int(src.width * scale)
        )

        band = src.read(
            1,
            out_shape=out_shape,
            resampling=Resampling.nearest
        )

        transform = src.transform * Affine.scale(
            src.width / out_shape[1],
            src.height / out_shape[0]
        )

        return band, transform, src.crs

def main_process(tile, datetime_range):

    # ============================================================
    # 1. STAC
    # ============================================================

    try:
        print(f"Iniciando Tile {tile}")

        catalog = pystac_client.Client.open(
            "https://data.inpe.br/bdc/stac/v1/"
        )

        search = catalog.search(
            collections=["S2-16D-2"],
            query={"bdc:tile": {"eq": tile}},
            datetime=datetime_range
        )

        item = next(search.items())


        # ============================================================
        # 2. LER BANDAS
        # ============================================================

        bands = ["B04", "B03", "B8A"]

        with ThreadPoolExecutor() as executor:

            results = list(
                executor.map(
                    read_band,
                    [item.assets[b].href for b in bands]
                )
            )

        red, transform, crs = results[0]
        green = results[1][0]
        nir = results[2][0]

        scene = np.stack([red, green, nir])


        # ============================================================
        # 3. OMNICLOUDMASK
        # ============================================================

        print(f"Omnicloud mask processing...")

        mask = omnicloudmask.predict_from_array(
            scene,
            inference_device="cpu",
            mosaic_device="cpu",
            inference_dtype="fp32",
            apply_no_data_mask=True,
            pred_classes=4,
        )

        mask_array = mask.squeeze()

        # ============================================================
        # 4. POLYGONIZE SOMENTE CLASSES 1, 2 E 3
        # ============================================================

        valid_mask = np.isin(
            mask_array,
            [1, 2, 3]
        )

        features = shapes(
            mask_array,
            mask=valid_mask,
            transform=transform
        )


        geoms = []
        classes = []

        for geom, value in features:

            geoms.append(shape(geom))
            classes.append(int(value))


        gdf = gpd.GeoDataFrame(
            {"class": classes},
            geometry=geoms,
            crs=crs
        )

        # ============================================================
        # 5. DISSOLVE GERAL DA OMNICLOUDMASK
        # ============================================================

        dissolved = unary_union(gdf.geometry)

        # Garantir somente Polygon
        if dissolved.geom_type == "Polygon":
            polygons = [dissolved]

        elif dissolved.geom_type == "MultiPolygon":
            polygons = list(dissolved.geoms)

        else:
            polygons = []

        gdf_dissolved = gpd.GeoDataFrame(
            geometry=polygons,
            crs=gdf.crs
        )

        result = gpd.GeoDataFrame(
            geometry=polygons,
            crs=gdf_dissolved.crs
        )

        # Remover vazias/nulas
        result = result[
            result.geometry.notna() &
            ~result.geometry.is_empty
        ].copy()

        result = result.to_crs("EPSG:4674")

        # ============================================================
        # 7. LOCALIZAR DIRETÓRIO DO TILE
        # ============================================================

        base_path = Path(
            "~/grupos/biomasbr-amazonia/sits-prodes/"
            "prodes.amz/data/class"
        ).expanduser()

        tile_path = base_path / tile
        post_processed = tile_path / "post_processed"


        processing_dir = next(
            p for p in post_processed.iterdir()
            if p.is_dir()
        )
        
        # ============================================================
        # 8. SALVAR GEOPACKAGE
        # ============================================================

        output_path = (
            processing_dir /
            f"rascunho_cloud_{tile}.gpkg"
        )

        result.to_file(
            output_path,
            layer="cloud_mask",
            driver="GPKG"
        )

        # Liberar memória
        try:
            del red, green, nir
            del scene, mask, mask_array, valid_mask
            del features, geoms, classes, gdf
            del dissolved, polygons, gdf_dissolved, result
        except NameError:
            pass
                
        print(f"Tile {tile} concluído!")
    
    finally:
        gc.collect()

#main_process("004015", datetime_range)

for tile in tiles:
     main_process(tile, datetime_range)