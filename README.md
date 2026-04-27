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
