# 📡 Muon Telescope Reconstruction & ML Framework

[![ROOT](https://img.shields.io/badge/ROOT-6.28-blue.svg)](https://root.cern/)
[![C++](https://img.shields.io/badge/C++-17-orange.svg)](https://isocpp.org/)
[![Thesis](https://img.shields.io/badge/Physics-Muon_Tomography-purple)]()

> **Context**  
> This repository houses the core analysis framework developed for my **Bachelor thesis in experimental particle physics**. 
The project tackles the reconstruction of atmospheric muon trajectories using a **Micromegas‑based tracking telescope**, with 
a strong emphasis on algorithmic noise rejection, calibration, and feature engineering for machine learning.

---

## 🔬 Telescope Hardware Configuration

The setup consists of **8 Micromegas detectors** arranged in a **4‑layer stereo configuration** (2 detectors per layer: `X` and `Y` 
coordinates), providing precise 2D spatial hit information.

| Layer | Height (Z) | Detector Pair |
|-------|------------|---------------|
| **1** (Bottom) | `0 mm` | `X1` & `Y1` |
| **2**           | `250 mm`   | `X2` & `Y2` |
| **3**           | `750 mm`   | `X3` & `Y3` |
| **4** (Top)    | `1000 mm`  | `X4` & `Y4` |

*(Normalized height `MM4y2 = 0.89` corresponds to the physical `1000 mm` height.)*

**Trigger System**  
Data acquisition is gated by two plastic scintillator paddles:
- **Bottom Scintillator**: `Z = 40 mm`
- **Top Scintillator**: `Z = 1130 mm`

---

## 🧩 Software Architecture & Dual Pipelines

The framework is split into **two distinct operational pipelines** to serve both traditional physics analysis and modern 
machine‑learning‑based classification. All code is written in **C++** using the **ROOT** data analysis framework.

### 📐 Pipeline 1: Standard Track Reconstruction (Physics Analysis)
This workflow generates physics plots (angular distributions, efficiency maps, residual analysis) from reconstructed tracks.

1. **`mmanalys.h`** – **Data Decoding**  
   Parses raw DAQ output, applies pedestal subtraction, and filters electronic noise specific to the Micromegas front‑end 
   electronics. Provides the `Loop()` method that iterates over events and fills hit collections.

2. **`histogram_booking.h`** – **Memory Management**  
   Books and maintains global pointers to all ROOT histograms and TTrees.

3. **`trackseeding.cpp`** – **Core Algorithm**  
   Implements a custom combinatorial seeding algorithm. This is the primary research contribution – it resolves multi‑hit ambiguities 
   and suppresses spurious hits caused by detector inefficiencies, enabling robust straight‑line fitting across the 4 layers.

4. **`plot.h`** – **Visualization**  
   Ingests the filled pointers to generate publication‑ready histograms.

---

### 🤖 Pipeline 2: Machine Learning Feature Extraction (ML Pipeline)
This pipeline reuses the same reconstruction logic but repurposes the output for supervised learning models. It is **executed 
in two passes** to first calibrate the spatial resolution empirically, then produce a clean feature set.

1. **`mmanalys.h`** (slightly modified) – **Data Decoding**  
   Used exactly as in Pipeline 1, but with a flag to switch between analytic and empirical error estimation.

2. **`histogram_booking.h`** – **Memory Management**  
   Shared with Pipeline 1.

3. **`micromegas.cpp`** – **Reconstruction + Feature Engineering** (two‑pass procedure):
   - **Pass 1** (`Loop(false)`) – processes events using analytic (geometric) error estimates. During this pass, the code fills 
     **unbiased residual histograms** (`h_res_vs_sizeX/Y`) vs cluster size. After the loop, it calls `FitSlicesY` and extracts 
     **empirical sigma values** for each cluster size (1–10 strips). These sigmas represent the true spatial resolution of the 
     detector for that cluster topology.
   - **Pass 2** (`Loop(true)`) – reprocesses the entire dataset, this time assigning per‑hit errors based on the empirical sigmas 
     from Pass 1. Now, for each event that forms a valid 4‑layer track, the code aggregates **feature vectors** (X, Y, Z positions, 
     total charge, cluster size per layer, and reconstructed track angles) and writes them to a structured **`ml_real_data.csv`** file. 
     This dataset is designed to train ML algorithms (e.g., BDTs, Neural Networks) to classify track quality or predict hit positions 
     with improved accuracy.

4. **`plot.h`** – **Diagnostic Visualization**  
   Used to validate the feature distributions (and the residual calibration) before feeding them into the ML training loop.

---

### 🔄 Data Flow Diagram
```mermaid
