# -*- coding: utf-8 -*-
"""
RiverMorph-DEM: 基于 DEM 的类 flowing-terrain 河流地貌生成网页工具

运行方式：
streamlit run river_morphology_app.py

功能：
1. 上传 DEM GeoTIFF；
2. 自动进行水文连通处理；
3. 计算 drainage / flow accumulation；
4. 提取河网；
5. 沿河道雕刻河谷；
6. 页面直接显示结果；
7. 下载 GeoTIFF、PNG 和 OBJ 三维模型。

适合：
- 河流地貌教学；
- DEM 驱动的河网与河谷生成；
- 后续论文方法原型。
"""

import io
import os
import zipfile
import tempfile
import heapq
from pathlib import Path

import numpy as np
import streamlit as st
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
from scipy import ndimage
import matplotlib.pyplot as plt
import trimesh
from PIL import Image


# =========================================================
# 页面设置
# =========================================================

st.set_page_config(
    page_title="RiverMorph-DEM",
    page_icon="🌊",
    layout="wide"
)


# =========================================================
# 基础函数
# =========================================================

def get_neighbors8(r, c, nrows, ncols):
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            rr = r + dr
            cc = c + dc
            if 0 <= rr < nrows and 0 <= cc < ncols:
                yield rr, cc


def fill_nodata_nearest(dem, nodata=None):
    arr = dem.copy().astype(np.float64)

    invalid = ~np.isfinite(arr)
    if nodata is not None:
        invalid |= (arr == nodata)

    if invalid.all():
        raise ValueError("DEM 全部为 NoData 或无效值，无法处理。")

    if invalid.any():
        indices = ndimage.distance_transform_edt(
            invalid,
            return_distances=False,
            return_indices=True
        )
        arr[invalid] = arr[tuple(indices[:, invalid])]

    return arr


def read_uploaded_dem(uploaded_file, max_size=800):
    """
    从 Streamlit 上传文件读取 DEM。
    """
    with tempfile.NamedTemporaryFile(delete=False, suffix=".tif") as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = tmp.name

    with rasterio.open(tmp_path) as src:
        height, width = src.height, src.width
        scale = max(height, width) / float(max_size)

        if scale > 1:
            out_height = int(round(height / scale))
            out_width = int(round(width / scale))
        else:
            out_height = height
            out_width = width

        dem = src.read(
            1,
            out_shape=(out_height, out_width),
            resampling=Resampling.bilinear
        ).astype(np.float64)

        transform = src.transform * Affine.scale(
            width / out_width,
            height / out_height
        )

        profile = src.profile.copy()
        profile.update(
            height=out_height,
            width=out_width,
            transform=transform,
            dtype="float32",
            count=1,
            compress="lzw"
        )

        nodata = src.nodata
        crs = src.crs

    os.remove(tmp_path)

    dem = fill_nodata_nearest(dem, nodata=nodata)

    return dem, profile, transform, crs, nodata


def priority_flood_receivers(dem):
    """
    Priority-Flood 风格的水文连通处理。

    从边界低处向内部推进，为每个像元记录下游接收像元 receiver。
    这一思想与 flowing-terrain 中“从低处向高处建立地形与排水关系”的逻辑相似。
    """
    nrows, ncols = dem.shape
    filled = dem.copy().astype(np.float64)

    visited = np.zeros((nrows, ncols), dtype=bool)
    receiver_r = np.full((nrows, ncols), -1, dtype=np.int32)
    receiver_c = np.full((nrows, ncols), -1, dtype=np.int32)

    heap = []
    counter = 0

    def push_cell(r, c, elev):
        nonlocal counter
        heapq.heappush(heap, (float(elev), counter, r, c))
        counter += 1

    # 边界作为出口
    for c in range(ncols):
        if not visited[0, c]:
            visited[0, c] = True
            push_cell(0, c, filled[0, c])
        if not visited[nrows - 1, c]:
            visited[nrows - 1, c] = True
            push_cell(nrows - 1, c, filled[nrows - 1, c])

    for r in range(nrows):
        if not visited[r, 0]:
            visited[r, 0] = True
            push_cell(r, 0, filled[r, 0])
        if not visited[r, ncols - 1]:
            visited[r, ncols - 1] = True
            push_cell(r, ncols - 1, filled[r, ncols - 1])

    visit_order = []

    while heap:
        elev, _, r, c = heapq.heappop(heap)
        visit_order.append((r, c))

        for rr, cc in get_neighbors8(r, c, nrows, ncols):
            if visited[rr, cc]:
                continue

            visited[rr, cc] = True

            if filled[rr, cc] < elev:
                filled[rr, cc] = elev

            receiver_r[rr, cc] = r
            receiver_c[rr, cc] = c

            push_cell(rr, cc, filled[rr, cc])

    return filled, receiver_r, receiver_c, visit_order


def compute_drainage_accumulation(receiver_r, receiver_c, visit_order):
    """
    从上游到下游累积 drainage / dampness。
    """
    nrows, ncols = receiver_r.shape
    acc = np.ones((nrows, ncols), dtype=np.float64)

    for r, c in reversed(visit_order):
        rr = receiver_r[r, c]
        cc = receiver_c[r, c]

        if rr >= 0 and cc >= 0:
            acc[rr, cc] += acc[r, c]

    return acc


def extract_streams(acc, percentile=99.3, min_size=8):
    threshold = np.percentile(acc, percentile)
    streams = acc >= threshold

    labels, num = ndimage.label(streams)
    if num > 0:
        counts = np.bincount(labels.ravel())
        keep = counts >= min_size
        keep[0] = False
        streams = keep[labels]

    return streams, threshold


def carve_valleys(
    dem,
    streams,
    transform,
    incision_depth=20.0,
    valley_width=150.0,
    smooth_sigma=0.6
):
    """
    沿河网进行高斯型下切。
    """
    cell_x = abs(transform.a)
    cell_y = abs(transform.e)
    cell_size = float(np.mean([cell_x, cell_y]))

    distance_pixels = ndimage.distance_transform_edt(~streams)
    distance_m = distance_pixels * cell_size

    sigma = max(valley_width, cell_size) / 2.0
    incision = incision_depth * np.exp(-(distance_m ** 2) / (2.0 * sigma ** 2))

    carved = dem - incision

    if smooth_sigma and smooth_sigma > 0:
        carved = ndimage.gaussian_filter(carved, sigma=smooth_sigma)

    return carved, incision, distance_m


def hillshade(dem, azimuth=315, altitude=45):
    x, y = np.gradient(dem)
    slope = np.pi / 2.0 - np.arctan(np.sqrt(x * x + y * y))
    aspect = np.arctan2(-x, y)

    az = np.deg2rad(azimuth)
    alt = np.deg2rad(altitude)

    shaded = (
        np.sin(alt) * np.sin(slope)
        + np.cos(alt) * np.cos(slope) * np.cos(az - aspect)
    )

    shaded = (shaded - shaded.min()) / (shaded.max() - shaded.min() + 1e-12)
    return shaded


def array_to_png_bytes(array, cmap="terrain", title=None, log_scale=False):
    arr = array.copy().astype(np.float64)

    if log_scale:
        arr = np.log1p(np.maximum(arr, 0))

    fig, ax = plt.subplots(figsize=(7, 6), dpi=180)
    im = ax.imshow(arr, cmap=cmap)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.036, pad=0.04)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def overlay_to_png_bytes(dem, acc, streams, title=None):
    hs = hillshade(dem)

    acc_log = np.log1p(acc)
    acc_norm = (acc_log - acc_log.min()) / (acc_log.max() - acc_log.min() + 1e-12)

    rgb = np.dstack([hs, hs, hs])

    # 红色表示 drainage potential
    rgb[..., 0] = np.maximum(rgb[..., 0], acc_norm)
    rgb[..., 1] = rgb[..., 1] * (1 - 0.45 * acc_norm)
    rgb[..., 2] = rgb[..., 2] * (1 - 0.45 * acc_norm)

    # 蓝色表示提取河网
    rgb[streams, 0] = 0.05
    rgb[streams, 1] = 0.35
    rgb[streams, 2] = 1.00

    fig, ax = plt.subplots(figsize=(8, 7), dpi=180)
    ax.imshow(np.clip(rgb, 0, 1))
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def save_geotiff_to_path(path, array, profile):
    profile_out = profile.copy()
    profile_out.update(dtype="float32", count=1, compress="lzw")

    with rasterio.open(path, "w", **profile_out) as dst:
        dst.write(array.astype(np.float32), 1)


def export_dem_to_obj_bytes(
    dem,
    transform,
    z_exaggeration=1.0,
    max_mesh_size=350
):
    """
    导出 OBJ，返回 bytes。
    """
    arr = dem.copy().astype(np.float64)
    nrows, ncols = arr.shape

    scale = max(nrows, ncols) / float(max_mesh_size)
    if scale > 1:
        out_rows = int(round(nrows / scale))
        out_cols = int(round(ncols / scale))
        zoom_r = out_rows / nrows
        zoom_c = out_cols / ncols
        arr = ndimage.zoom(arr, (zoom_r, zoom_c), order=1)

        transform = transform * Affine.scale(ncols / out_cols, nrows / out_rows)
        nrows, ncols = arr.shape

    cell_x = abs(transform.a)
    cell_y = abs(transform.e)

    x = np.arange(ncols) * cell_x
    y = np.arange(nrows) * cell_y
    xx, yy = np.meshgrid(x, y)

    z = arr * z_exaggeration

    vertices = np.column_stack([
        xx.ravel(),
        yy.ravel(),
        z.ravel()
    ])

    faces = []
    for r in range(nrows - 1):
        for c in range(ncols - 1):
            v0 = r * ncols + c
            v1 = r * ncols + c + 1
            v2 = (r + 1) * ncols + c
            v3 = (r + 1) * ncols + c + 1

            faces.append([v0, v2, v1])
            faces.append([v1, v2, v3])

    mesh = trimesh.Trimesh(
        vertices=vertices,
        faces=np.asarray(faces),
        process=False
    )

    with tempfile.NamedTemporaryFile(delete=False, suffix=".obj") as tmp:
        obj_path = tmp.name

    mesh.export(obj_path)

    with open(obj_path, "rb") as f:
        obj_bytes = f.read()

    os.remove(obj_path)
    return obj_bytes


def build_result_zip(
    profile,
    dem,
    filled,
    acc,
    streams,
    carved,
    incision,
    obj_bytes,
    overlay_png,
    carved_overlay_png
):
    """
    把主要结果打包成 zip。
    """
    tmp_dir = tempfile.mkdtemp()

    paths = {}

    paths["01_input_dem.tif"] = os.path.join(tmp_dir, "01_input_dem.tif")
    paths["02_filled_dem.tif"] = os.path.join(tmp_dir, "02_filled_dem.tif")
    paths["03_drainage_accumulation.tif"] = os.path.join(tmp_dir, "03_drainage_accumulation.tif")
    paths["04_stream_network.tif"] = os.path.join(tmp_dir, "04_stream_network.tif")
    paths["05_carved_dem.tif"] = os.path.join(tmp_dir, "05_carved_dem.tif")
    paths["06_incision_depth.tif"] = os.path.join(tmp_dir, "06_incision_depth.tif")
    paths["07_carved_terrain.obj"] = os.path.join(tmp_dir, "07_carved_terrain.obj")
    paths["08_flowing_terrain_like_overlay.png"] = os.path.join(tmp_dir, "08_flowing_terrain_like_overlay.png")
    paths["09_carved_overlay.png"] = os.path.join(tmp_dir, "09_carved_overlay.png")

    save_geotiff_to_path(paths["01_input_dem.tif"], dem, profile)
    save_geotiff_to_path(paths["02_filled_dem.tif"], filled, profile)
    save_geotiff_to_path(paths["03_drainage_accumulation.tif"], acc, profile)
    save_geotiff_to_path(paths["04_stream_network.tif"], streams.astype(np.float32), profile)
    save_geotiff_to_path(paths["05_carved_dem.tif"], carved, profile)
    save_geotiff_to_path(paths["06_incision_depth.tif"], incision, profile)

    with open(paths["07_carved_terrain.obj"], "wb") as f:
        f.write(obj_bytes)

    with open(paths["08_flowing_terrain_like_overlay.png"], "wb") as f:
        f.write(overlay_png.getvalue())

    with open(paths["09_carved_overlay.png"], "wb") as f:
        f.write(carved_overlay_png.getvalue())

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, path in paths.items():
            zf.write(path, arcname=name)

    zip_buffer.seek(0)
    return zip_buffer


# =========================================================
# 页面主体
# =========================================================

st.title("🌊 RiverMorph-DEM：基于 DEM 的河流地貌生成工具")

st.markdown(
    """
本工具用于从输入 DEM 中生成类似 **flowing-terrain** 的河流汇流潜势、河网和三维河谷地形。
核心思想是：从真实 DEM 出发，建立水文连通关系，计算 drainage accumulation，并沿高汇流区雕刻河谷。
"""
)

with st.sidebar:
    st.header("参数设置")

    uploaded_dem = st.file_uploader(
        "上传 DEM GeoTIFF 文件",
        type=["tif", "tiff"]
    )

    max_size = st.slider(
        "参与计算的最大 DEM 尺寸",
        min_value=300,
        max_value=1500,
        value=800,
        step=100,
        help="数值越大结果越精细，但计算越慢。"
    )

    stream_percentile = st.slider(
        "河网提取阈值百分位",
        min_value=95.0,
        max_value=99.9,
        value=99.3,
        step=0.1,
        help="数值越高，河网越稀疏；数值越低，支流越多。"
    )

    min_stream_size = st.slider(
        "最小河网斑块大小",
        min_value=1,
        max_value=100,
        value=8,
        step=1
    )

    incision_depth = st.slider(
        "最大河谷下切深度",
        min_value=1.0,
        max_value=100.0,
        value=20.0,
        step=1.0,
        help="单位与 DEM 高程单位一致，通常为米。"
    )

    valley_width = st.slider(
        "河谷宽度控制参数",
        min_value=10.0,
        max_value=1000.0,
        value=150.0,
        step=10.0,
        help="单位与 DEM 水平单位一致，投影坐标下通常为米。"
    )

    smooth_sigma = st.slider(
        "地形平滑强度",
        min_value=0.0,
        max_value=5.0,
        value=0.6,
        step=0.1
    )

    z_exaggeration = st.slider(
        "三维模型高程夸张系数",
        min_value=0.2,
        max_value=5.0,
        value=1.5,
        step=0.1
    )

    max_mesh_size = st.slider(
        "OBJ 模型最大网格尺寸",
        min_value=100,
        max_value=800,
        value=350,
        step=50,
        help="越大模型越精细，但 OBJ 文件越大。"
    )

    run_button = st.button("开始生成", type="primary")


if uploaded_dem is None:
    st.info("请在左侧上传一个 DEM GeoTIFF 文件，然后点击“开始生成”。")
    st.stop()


if run_button:
    try:
        progress = st.progress(0)
        status = st.empty()

        status.write("正在读取 DEM...")
        dem, profile, transform, crs, nodata = read_uploaded_dem(
            uploaded_dem,
            max_size=max_size
        )
        progress.progress(10)

        status.write("正在进行水文连通处理...")
        filled, receiver_r, receiver_c, visit_order = priority_flood_receivers(dem)
        progress.progress(30)

        status.write("正在计算 drainage / flow accumulation...")
        acc = compute_drainage_accumulation(receiver_r, receiver_c, visit_order)
        progress.progress(45)

        status.write("正在提取河网...")
        streams, threshold = extract_streams(
            acc,
            percentile=stream_percentile,
            min_size=min_stream_size
        )
        progress.progress(60)

        status.write("正在雕刻河谷...")
        carved, incision, distance_m = carve_valleys(
            dem=filled,
            streams=streams,
            transform=transform,
            incision_depth=incision_depth,
            valley_width=valley_width,
            smooth_sigma=smooth_sigma
        )
        progress.progress(75)

        status.write("正在生成图像和三维模型...")
        input_dem_png = array_to_png_bytes(dem, cmap="terrain", title="Input DEM")
        acc_png = array_to_png_bytes(acc, cmap="inferno", title="Drainage accumulation", log_scale=True)
        stream_png = array_to_png_bytes(streams.astype(float), cmap="Blues", title="Extracted stream network")
        carved_png = array_to_png_bytes(carved, cmap="terrain", title="Carved DEM")
        incision_png = array_to_png_bytes(incision, cmap="magma", title="Incision depth")

        overlay_png = overlay_to_png_bytes(
            filled,
            acc,
            streams,
            title="Flowing-terrain-like drainage map"
        )

        carved_overlay_png = overlay_to_png_bytes(
            carved,
            acc,
            streams,
            title="Carved terrain + drainage + streams"
        )

        obj_bytes = export_dem_to_obj_bytes(
            carved,
            transform,
            z_exaggeration=z_exaggeration,
            max_mesh_size=max_mesh_size
        )

        zip_buffer = build_result_zip(
            profile=profile,
            dem=dem,
            filled=filled,
            acc=acc,
            streams=streams,
            carved=carved,
            incision=incision,
            obj_bytes=obj_bytes,
            overlay_png=overlay_png,
            carved_overlay_png=carved_overlay_png
        )

        progress.progress(100)
        status.success("生成完成！")

        st.subheader("基本信息")

        col_info1, col_info2, col_info3, col_info4 = st.columns(4)
        col_info1.metric("DEM 行数", dem.shape[0])
        col_info2.metric("DEM 列数", dem.shape[1])
        col_info3.metric("河网阈值", f"{threshold:.2f}")
        col_info4.metric("河网像元数", int(streams.sum()))

        st.subheader("结果预览")

        tab1, tab2, tab3, tab4, tab5 = st.tabs(
            [
                "原始 DEM",
                "汇流累积",
                "河网",
                "类 flowing-terrain 效果",
                "雕刻后地形"
            ]
        )

        with tab1:
            st.image(input_dem_png, caption="输入 DEM", use_container_width=True)

        with tab2:
            st.image(acc_png, caption="Drainage accumulation，对数尺度", use_container_width=True)

        with tab3:
            st.image(stream_png, caption="提取出的河网", use_container_width=True)

        with tab4:
            st.image(
                overlay_png,
                caption="灰度为地形阴影，红色为汇流潜势，蓝色为提取河网",
                use_container_width=True
            )

        with tab5:
            col_a, col_b = st.columns(2)
            with col_a:
                st.image(carved_png, caption="河谷雕刻后的 DEM", use_container_width=True)
            with col_b:
                st.image(incision_png, caption="河谷下切深度", use_container_width=True)

            st.image(
                carved_overlay_png,
                caption="雕刻后地形 + 汇流潜势 + 河网",
                use_container_width=True
            )

        st.subheader("结果下载")

        col_d1, col_d2, col_d3 = st.columns(3)

        with col_d1:
            st.download_button(
                label="下载全部结果 ZIP",
                data=zip_buffer,
                file_name="RiverMorph_DEM_results.zip",
                mime="application/zip"
            )

        with col_d2:
            st.download_button(
                label="下载三维 OBJ 模型",
                data=obj_bytes,
                file_name="carved_terrain.obj",
                mime="text/plain"
            )

        with col_d3:
            st.download_button(
                label="下载效果图 PNG",
                data=carved_overlay_png,
                file_name="carved_overlay.png",
                mime="image/png"
            )

        st.subheader("参数解释")

        st.markdown(
            f"""
本次运行中，河网提取采用汇流累积的 **{stream_percentile:.1f}% 百分位阈值**。
阈值越高，保留的河道越少，通常更接近主干河；阈值越低，支流数量越多。

河谷雕刻采用高斯型下切函数。最大下切深度为 **{incision_depth:.1f}**，
河谷宽度控制参数为 **{valley_width:.1f}**。因此，距离河道越近的像元下切越强，
距离河道越远的像元影响逐渐减弱。
"""
        )

    except Exception as e:
        st.error(f"运行出错：{e}")
        st.exception(e)

else:
    st.info("DEM 已上传。请在左侧调整参数，然后点击“开始生成”。")
