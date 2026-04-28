# RiverMorphDEM

**RiverMorphDEM** is a DEM-constrained river morphology generation tool. It uses a Digital Elevation Model (DEM) as input and automatically generates drainage accumulation, stream networks, carved river valleys, and exportable three-dimensional terrain models.

The tool is designed for river geomorphology teaching, terrain cognition, procedural landscape generation, and preliminary research on DEM-based river morphology modeling.

中文简介：  
**RiverMorphDEM** 是一个基于数字高程模型（DEM）的河流地貌生成工具，可从输入 DEM 中自动计算汇流潜势、提取河网、雕刻河谷，并导出可视化图像、GeoTIFF 栅格结果和 OBJ 三维实体模型。该工具适用于河流地貌教学、课程作业、地形认知训练和程序化地貌建模研究。

---

## 1. Main Features

RiverMorphDEM currently supports the following functions:

- Upload GeoTIFF DEM data through a web interface.
- Read and preprocess DEM files.
- Perform Priority-Flood-style hydrological conditioning.
- Compute drainage accumulation.
- Extract stream networks based on drainage thresholds.
- Generate flowing-terrain-like drainage visualization.
- Carve river valleys along extracted stream networks.
- Preview carved terrain in an interactive 3D Plotly viewer (downsampled for speed).
- Export processed DEM, stream network, drainage accumulation, incision depth, PNG figures, and OBJ terrain models.
- Download all results as a ZIP file.

---

## 2. Method Overview

The workflow of RiverMorphDEM is as follows:

```text
Input DEM
   ↓
NoData handling and optional resampling
   ↓
Hydrological conditioning
   ↓
Drainage accumulation calculation
   ↓
Stream network extraction
   ↓
River valley carving
   ↓
2D visualization and 3D mesh generation
   ↓
GeoTIFF / PNG / OBJ / ZIP output
```

---

## 3. Environment Requirements

Install dependencies from `requirements.txt`:

```bash
pip install -r requirements.txt
```

The project currently depends on:

- streamlit
- numpy
- scipy
- rasterio
- matplotlib
- trimesh
- pillow
- plotly

---

## 4. Run the App Locally

From the repository root, start Streamlit with:

```bash
python -m streamlit run app.py
```

After startup, open the local URL shown in the terminal (usually `http://localhost:8501`).

---

## 5. Quick Start / 快速开始

```bash
# 1) 克隆仓库
git clone <your-repo-url>
cd RiverMorphDEM

# 2) 安装依赖
pip install -r requirements.txt

# 3) 运行 Streamlit 应用
python -m streamlit run app.py
```

Quick workflow in browser:

1. Upload a DEM GeoTIFF file.
2. Select a preset in **Parameter preset / 参数预设** (or keep custom).
3. Optionally tune sliders for stream extraction and valley carving.
4. Click **开始生成**.
5. Review 2D + 3D previews and download ZIP/OBJ outputs.

---

## 6. Limitations / 局限性

- 当前工具定位为教学与研究原型（teaching/research prototype）。
- 河网提取与河谷雕刻效果依赖 DEM 数据质量和参数设置。
- 不建议直接用于工程设计、防洪评价或高精度水动力模拟。
- 为保证网页运行效率，较大 DEM 会在计算或预览阶段进行降采样处理。

---

## 7. Input Data Requirements & Common Problems

After uploading a DEM, the app now provides a **DEM quality diagnostics** panel, including:

- Raster rows and columns.
- Elevation minimum and maximum.
- NoData value and NoData pixel ratio.
- Pixel resolution (X/Y).
- CRS string and whether it is projected.
- A quick judgement on whether large invalid regions may exist.

Quality warnings shown in the page:

1. If the DEM is in geographic coordinates (longitude/latitude), the app warns that valley width should not be directly interpreted as meters and recommends reprojecting to a metric CRS first.
2. If NoData ratio is high (for example >10%), the app warns to crop valid terrain first to reduce interpolation artifacts near boundaries.
3. If elevation range appears abnormal (too small or with extreme outliers), the app shows an elevation-range hint so users can inspect data quality before processing.

---

## 8. Parameter Presets / 参数预设

The sidebar includes a preset selector: **Parameter preset / 参数预设**.

Available presets and typical use cases:

- **UAV / local gully DEM**: for UAV-scale data, local gullies, and short channel segments.
- **Small catchment DEM**: for small watershed or hillslope-channel systems.
- **Large basin demo**: for large basin classroom demos or regional DEM examples.

Preset values:

| Preset | max_size | stream_percentile | min_stream_size | incision_depth | valley_width | smooth_sigma | z_exaggeration | max_mesh_size |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| UAV / local gully DEM | 1200 | 99.5 | 30 | 3 | 8 | 0.3 | 1.5 | 500 |
| Small catchment DEM | 1000 | 99.3 | 20 | 10 | 50 | 0.5 | 1.5 | 450 |
| Large basin demo | 800 | 99.0 | 10 | 20 | 150 | 0.6 | 1.5 | 350 |

You can still manually adjust all sliders after selecting a preset.

---

## 9. 3D Preview / 三维预览

After result generation, the app shows an interactive **3D Preview / 三维预览** based on the **carved DEM**:

- Implemented with Plotly (`go.Surface`) in-page preview.
- Supports rotation and zoom for quick terrain-shape inspection.
- Uses elevation as the Z value and applies the current `z_exaggeration`.
- Uses downsampling (preview max size around 200–300 cells) to keep the page responsive.

> Note: the 3D preview is only for quick visualization. For full-resolution geometry, download the OBJ output.
