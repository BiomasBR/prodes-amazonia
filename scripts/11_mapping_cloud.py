from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import gc
import geopandas as gpd
import numpy as np
import omnicloudmask
import pystac_client
import rasterio
import torch
import zipfile

from affine import Affine
from rasterio.enums import Resampling
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union

import os
import re
import sys
import pyproj

# Todos os tiles abaixo serão processados automaticamente.
TILES = [
'002015',
'005012',
'005013',
'006008',
'006010',
'006011',
'007003',
'007004',
'007005',
'007011',
'008004',
'008012',
'009005',
'009006',
'009011',
'011001',
'011007',
'011008',
'011009',
'012006',
'013007',
'014008',
'015000',
'031008'
]

# Data de cada tipo de nuvem, valendo para TODOS os tiles.
# Aceita "AAAA-MM-DD" ou intervalo "AAAA-MM-DD/AAAA-MM-DD".
DATAS = {
    "primaria":   "2026-07-28/2026-07-30",             
    "secundaria": "2026-08-13/2026-08-15",
}

# Exceções opcionais por tile (sobrescrevem DATAS só nas chaves informadas).
# Exemplo: {"017006": {"secundaria": "2026-07-20/2026-07-23"}}
DATAS_POR_TILE = {}

OUTPUT_DIR = Path(
    "~/grupos/biomasbr-amazonia/sits-prodes/"
            "prodes.amz/data/raw/cloud_np"
).expanduser()

EXECUTAR_GERACAO = True            
EXECUTAR_POS = True                
PULAR_GPKG_EXISTENTES = False    

TORCH_THREADS = 16
CRS_FINAL = "EPSG:4674"

# Nomes dos arquivos de cada tipo de nuvem (um .zip separado por tipo e por tile)
TIPOS = {
    "primaria": {
        "gpkg": "rascunho_cloud_{tile}.gpkg",
        "zip": "rascunho_cloud_{tile}.zip",
    },
    "secundaria": {
        "gpkg": "rascunho_cloud_secundaria_{tile}.gpkg",
        "zip": "rascunho_cloud_secundaria_{tile}.zip",
    },
}

EXTENSOES = ["dbf", "prj", "shp", "shx"]        # as mesmas do script em R

# ============================================================
# CONFIGURAÇÃO PROJ
# ============================================================

conda_proj_path = "/opt/conda/envs/geospatial/share/proj"

if os.path.exists(conda_proj_path):
    os.environ["PROJ_DATA"] = conda_proj_path
    pyproj.datadir.set_data_dir(conda_proj_path)
else:
    pyproj_datadir = pyproj.datadir.get_data_dir()
    os.environ["PROJ_DATA"] = pyproj_datadir
    pyproj.datadir.set_data_dir(pyproj_datadir)

torch.set_num_threads(TORCH_THREADS)

# ============================================================
# VALIDAÇÃO DA CONFIGURAÇÃO
# ============================================================

def normalizar_datas(valor, rotulo):

    if isinstance(valor, (list, tuple)):
        valor = "/".join(str(v) for v in valor)

    partes = [p.strip() for p in str(valor).split("/")]

    if len(partes) == 1:
        partes = partes * 2

    if len(partes) != 2:
        raise ValueError(f"{rotulo}: data inválida {valor!r}")

    try:
        ini, fim = (datetime.strptime(p, "%Y-%m-%d") for p in partes)
    except ValueError:
        raise ValueError(
            f"{rotulo}: data inválida {valor!r} "
            f"(use AAAA-MM-DD ou AAAA-MM-DD/AAAA-MM-DD)"
        ) from None

    if ini > fim:
        raise ValueError(f"{rotulo}: data inicial maior que a final ({valor!r})")

    return f"{ini:%Y-%m-%d}/{fim:%Y-%m-%d}"

def montar_jobs():
    """Valida toda a configuração antes de qualquer processamento pesado.
    Retorna {tile: {tipo: datetime_range}}."""

    tiles = list(dict.fromkeys(str(t).strip() for t in TILES))

    erros = [f"tile inválido: {t!r}" for t in tiles if not re.fullmatch(r"\d{6}", t)]
    erros += [
        f"DATAS_POR_TILE contém tile fora de TILES: {t}"
        for t in DATAS_POR_TILE if t not in tiles
    ]

    jobs = {}

    for tile in tiles:

        datas = {**DATAS, **DATAS_POR_TILE.get(tile, {})}

        for tipo in TIPOS:

            valor = datas.get(tipo)

            if valor is None:
                continue

            try:
                jobs.setdefault(tile, {})[tipo] = normalizar_datas(
                    valor, f"{tile}/{tipo}"
                )
            except ValueError as e:
                erros.append(str(e))

    if erros:
        raise SystemExit("Configuração inválida:\n  - " + "\n  - ".join(erros))

    if not jobs:
        raise SystemExit("Nenhum job configurado (DATAS vazio).")

    return jobs

# ============================================================
# LEITURA DAS BANDAS
# ============================================================

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

# ============================================================
# ETAPA 1 - PROCESSAMENTO (GERAÇÃO DA NUVEM)
# ============================================================

def main_process(tile, datetime_range, output_path):

    try:

        print(f"Iniciando Tile {tile}")

        # ====================================================
        # 1. STAC
        # ====================================================

        catalog = pystac_client.Client.open(
            "https://data.inpe.br/bdc/stac/v1/"
        )

        search = catalog.search(
            collections=["S2-16D-2"],
            query={"bdc:tile": {"eq": tile}},
            datetime=datetime_range
        )

        item = next(search.items(), None)

        if item is None:
            raise LookupError(
                f"nenhuma cena em S2-16D-2 para o tile {tile} em {datetime_range}"
            )

        # ====================================================
        # 2. LER BANDAS
        # ====================================================

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

        scene = np.stack([
            red,
            green,
            nir
        ])

        # ====================================================
        # 3. OMNICLOUDMASK
        # ====================================================

        print(f"Omnicloud mask processing - Tile {tile}...")

        mask = omnicloudmask.predict_from_array(
            scene,
            inference_device="cpu",
            mosaic_device="cpu",
            inference_dtype="fp32",
            apply_no_data_mask=True,
            pred_classes=4,
        )

        mask_array = mask.squeeze()

        # ====================================================
        # 4. POLYGONIZE SOMENTE CLASSES 1, 2 E 3
        # ====================================================

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

        # ====================================================
        # 5. DISSOLVE GERAL DA OMNICLOUDMASK
        # ====================================================

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

        # ====================================================
        # 6. REMOVER GEOMETRIAS VAZIAS / NULAS
        # ====================================================

        result = result[
            result.geometry.notna() &
            ~result.geometry.is_empty
        ].copy()

        # ====================================================
        # 7. DIRETÓRIO ÚNICO DE SAÍDA
        # ====================================================

        output_dir = OUTPUT_DIR

        output_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        # ====================================================
        # 8. SALVAR GEOPACKAGE
        # ====================================================

        if output_path.exists():
            output_path.unlink()

        result.to_file(
            output_path,
            layer="cloud_mask",
            driver="GPKG",
            engine="pyogrio",
            mode="w",
            geometry_type="Polygon"
        )

        print(
            f"Tile {tile} concluído! "
            f"Resultado salvo em: {output_path}"
        )

        # ====================================================
        # 9. LIBERAR MEMÓRIA
        # ====================================================

        try:

            del red
            del green
            del nir
            del scene
            del mask
            del mask_array
            del valid_mask
            del features
            del geoms
            del classes
            del gdf
            del dissolved
            del polygons
            del gdf_dissolved
            del result

        except NameError:
            pass


    finally:

        gc.collect()

# ============================================================
# ETAPA 2 - PÓS-PROCESSAMENTO
# ============================================================

def caminho_gpkg(tile, tipo):
    return OUTPUT_DIR / TIPOS[tipo]["gpkg"].format(tile=tile)

def remover_shapefile(shp_path):
    """Equivalente ao delete_dsn = TRUE do st_write: apaga o shapefile e seus arquivos."""
    for ext in EXTENSOES + ["cpg"]:
        shp_path.with_suffix("." + ext).unlink(missing_ok=True)

def gpkg_para_shapefile(gpkg_path):

    shp_path = gpkg_path.with_suffix(".shp")

    file = gpd.read_file(gpkg_path, engine="pyogrio")

    # st_cast("POLYGON")
    file = file.explode(index_parts=False, ignore_index=True)

    # file$id <- as.integer(seq_len(nrow(file)))
    file["id"] = np.arange(1, len(file) + 1, dtype="int32")

    # select("id") |> st_transform("EPSG:4674")
    file = file[["id", "geometry"]].to_crs(CRS_FINAL)

    remover_shapefile(shp_path)

    file.to_file(
        shp_path,
        driver="ESRI Shapefile",
        engine="pyogrio",
        geometry_type="Polygon"
    )

    return shp_path

def pos_processar_tipo(tile, tipo):
    """gpkg -> shapefile -> zip (arquivos na raiz do zip) -> remove o gpkg compactado."""

    gpkg = caminho_gpkg(tile, tipo)

    if not gpkg.exists():
        raise FileNotFoundError(f"falta {gpkg.name}")

    zip_path = OUTPUT_DIR / TIPOS[tipo]["zip"].format(tile=tile)
    shp_path = None

    try:

        # 1) conversão do .gpkg
        shp_path = gpkg_para_shapefile(gpkg)

        arquivos = [shp_path.with_suffix("." + ext) for ext in EXTENSOES]

        for arq in arquivos:
            if not arq.exists():
                raise FileNotFoundError(f"shapefile incompleto: falta {arq.name}")

        # 2) zip (remove o anterior, se existir)
        if zip_path.exists():
            zip_path.unlink()

        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for arq in arquivos:
                zf.write(arq, arq.name)

        # 3) verificação de integridade antes de apagar qualquer coisa
        with zipfile.ZipFile(zip_path) as zf:
            if zf.testzip() is not None or set(zf.namelist()) != {a.name for a in arquivos}:
                raise RuntimeError("zip gerado não passou na verificação")

        # 4) zip validado: remove o .gpkg que foi compactado
        gpkg.unlink(missing_ok=True)

    except Exception:

        # não deixa zip parcial; o .gpkg continua na pasta
        zip_path.unlink(missing_ok=True)
        raise

    finally:

        # shapefile intermediário nunca fica na pasta
        if shp_path is not None:
            remover_shapefile(shp_path)

        gc.collect()

    return zip_path

# ============================================================
# EXECUÇÃO
# ============================================================

def main():

    jobs = montar_jobs()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    n_tiles = len(jobs)
    resumo = {}

    for i, (tile, datas) in enumerate(jobs.items(), start=1):

        print(f"\n===== Tile {tile} ({i}/{n_tiles}) =====")

        status = {tipo: "-" for tipo in TIPOS}
        pos = {tipo: "-" for tipo in TIPOS}

        # ---------------- Etapa 1: geração ----------------
        if EXECUTAR_GERACAO:

            for tipo, datetime_range in datas.items():

                gpkg = caminho_gpkg(tile, tipo)

                if PULAR_GPKG_EXISTENTES and gpkg.exists():
                    status[tipo] = "existente"
                    print(f"{gpkg.name} já existe - pulando")
                    continue

                print(f"-> Nuvem {tipo} | tile {tile} | {datetime_range}")

                try:
                    main_process(tile, datetime_range, gpkg)
                    status[tipo] = "ok"

                except Exception as e:
                    status[tipo] = "erro"
                    print(f"ERROR in tile: {tile} ({tipo}): {e}")

        # ---------------- Etapas 2 e 3: pós + zip + limpeza (por tipo) ----------------
        if EXECUTAR_POS:

            for tipo in datas:

                if status[tipo] == "erro":
                    pos[tipo] = "geração falhou"
                    continue

                try:
                    zip_path = pos_processar_tipo(tile, tipo)
                    pos[tipo] = "ok"
                    print(f"{zip_path.name} criado; {caminho_gpkg(tile, tipo).name} removido")

                except Exception as e:
                    pos[tipo] = "erro"
                    print(f"ERROR in tile: {tile} ({tipo}, pós-processamento): {e}")

        resumo[tile] = (status, pos)

    # ---------------- Resumo ----------------
    print("\n" + "=" * 78)
    print(f"{'TILE':<8}{'GERAÇÃO prim.':<15}{'GERAÇÃO sec.':<15}{'ZIP prim.':<16}{'ZIP sec.'}")
    print("-" * 78)

    falhas = 0

    for tile, (status, pos) in resumo.items():
        print(
            f"{tile:<8}{status['primaria']:<15}{status['secundaria']:<15}"
            f"{pos['primaria']:<16}{pos['secundaria']}"
        )
        falhas += sum(v in ("erro", "geração falhou") for v in [*status.values(), *pos.values()]) 

    print("=" * 78)
    print(f"Tiles: {n_tiles} | falhas: {falhas}")

    return 1 if falhas else 0

if __name__ == "__main__":
    sys.exit(main())
