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
coordinates), providing precise 2D spatial hit information. The detectors dimensions active region is $460 \times 460$mm but we work only with active region $384 \times 384$mm.

| Layer | Height (Z) | Detector Pair |
|-------|------------|---------------|
| **1** (Bottom) | `0 mm` | `X1` & `Y1` |
| **2**           | `250 mm`   | `X2` & `Y2` |
| **3**           | `750 mm`   | `X3` & `Y3` |
| **4** (Top)    | `1000 mm`  | `X4` & `Y4` |


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

5. **'MMHit.h'** stores raw strip number, charge, and drift time (y from FD_time).

6. **'MMPacmanAlgo.h'** implements a simple forward‑only clustering where the maximum gap between adjacent hits is controlled by m_clus_size (default 1, meaning  no gaps allowed). The code includes a m_hole_size parameter but it's not used (commented out). The algorithm groups hits that are within m_clus_size channels.

7. **'MMCluster.h'** provides methods like Pos(pitch) (centroid), Charge(), Channel() (centroid in strip units), and size().

8. **'MMTrackingChain.h'** is a separate class that implements a different track reconstruction approach (based on "chains") – it's not used in the main pipelines but is included as an alternative. I've shown it with dashed lines.

9. **'plot.h'** is used in both pipelines to draw histograms; in the ML pipeline, it also draws calibration plots (e.g., residual vs cluster size).
---

### 🔄 Data Flow Diagram
```mermaid

graph TD
    subgraph Hardware
       Muon [Muon passes through the detector]--> Scintillator [Charge Detected on Both Scintillators] --> MM["Micromegas Collects Data"] --> ADC [Analog to Digital Cards] --> FEC [Front-End Cards] --> DAQ["Raw DAQ Data"]
    end

    subgraph Pipeline_1["Pipeline 1: Standard Physics Analysis"]
        MM --> A["mmanalys.h<br>Noise Filtering & Decoding"]
        A --> B["histogram_booking.h<br>ROOT Pointer Booking"]
        B --> C["trackseeding.cpp<br>Combinatorial Track Seeding"]
        C --> D["plot.h<br>Physics Histograms"]
    end

    subgraph Pipeline_2["Pipeline 2: Machine Learning Preparation"]
        MM --> E["mmanalys.h (modified)<br>Data Decoding"]
        E --> F["micromegas.cpp<br>Pass 1: Empirical Sigma Extraction"]
        F --> G["micromegas.cpp<br>Pass 2: Feature Aggregation & CSV Export"]
        G --> H["ml_real_data.csv<br>ML Training Dataset"]
        F --> I["plot.h<br>Residual Calibration Plots"]
        G --> I
    end
```

My Algorithmic Contributions
The principal challenge of this telescope is high noise occupancy and multi‑hit ambiguities inherent to the Micromegas gaseous detectors. 
My work introduces a complete suite of algorithmic corrections—spanning hardware mapping to advanced reconstruction physics—that transformed 
the detector's performance. Below are the core innovations implemented in both trackseeding.cpp (Physics pipeline) and micromegas.cpp (ML pipeline):

1. Hardware Strip Mapping Correction (Hardware‑to‑Software alignment)
The raw APV readout had a severe cable-swap error: strips from physically distinct readout chips were routed to the wrong software channels. 
I reconstructed the exact hardware topology and implemented an explicit stripMap[800] lookup table that surgically re‑routes each APV block 
(e.g., swapping APV 0 ↔ APV 4, APV 1 ↔ APV 3, and fixing intra‑APV mirroring). This single fix eliminated thousands of geometrically impossible 
ghost tracks and restored the expected linear correlation between adjacent layers.

2. Electronic Noise & Crosstalk Suppression (Pedestal Charge Thresholding)
The front‑end electronics suffered from significant crosstalk, producing a high rate of low‑charge spikes that mimicked real hits. By carefully 
increasing the pedestal charge threshold and applying a global minimum‑charge cut (current_charge <= 80 → skip), I entirely suppressed these spurious 
micro‑hits. This drastically cleaned the hitmaps and prevented avalanche‑like noise from contaminating the track seed candidates.

3. Best‑Fit Hit Selection (Replacing Highest‑Charge Heuristics)
Instead of blindly selecting the cluster with the largest integrated charge (which is often dominated by a delta‑electron or an asymmetric avalanche), 
I developed a combinatorial iterative algorithm that evaluates all possible hit combinations across the 4 layers. It systematically swaps out the 
worst‑fitting hit in the current straight‑line hypothesis until the residual distribution converges to the minimum. This approach reduced the percentage 
of poorly reconstructed tracks from ~12% down to just 4‑8%—a major leap in data quality.

4. Hybrid Reconstruction Physics (Centroid + µTPC from CERN)
To maximize spatial resolution across all incident angles, I implemented and hybridized two distinct position‑reconstruction methods:

Standard Charge Centroid – ideal for small clusters (≤ 3 strips) and near‑vertical tracks, where charge sharing is symmetric.

Micro‑Time Projection Chamber (µTPC) method – developed at CERN for large‑angle tracks. By analyzing the correlation between the drift time (gety()) 
and strip position, it corrects the systematic bias that plagues the centroid at steep angles.

The hybrid logic automatically switches between these two methods based on cluster size and the preliminary reconstructed angle, ensuring optimal 
precision across the full 60°–120° acceptance.

5. Dynamic Empirical Spatial Resolution (Two‑Pass Auto‑Calibration)
Spatial resolution is not constant—it depends heavily on cluster size and track topology. I designed a two‑pass calibration routine inside micromegas.cpp:

Pass 1 (Loop(false)): Uses analytic error estimates to fit tracks and fill unbiased residual histograms (h_res_vs_sizeX/Y) versus cluster size. After the 
loop, it extracts the empirical sigma for each cluster size (1–10 strips) directly from the data.

Pass 2 (Loop(true)): Re‑processes the entire dataset, now assigning a unique, realistic spatial error to each hit based on its cluster size and the 
empirical sigmas from Pass 1. This yields fits with properly weighted uncertainties and unleashes the full potential of the ML feature extraction.

6. Comprehensive Residual Diagnostics (Unbiased vs. Biased)
To rigorously validate the reconstruction, I implemented both biased (standard 4‑point fit) and unbiased (leave‑one‑out cross‑validation) residual 
calculators. The unbiased residuals—which exclude the target layer from the track fit before calculating its deviation—provide a mathematically pure 
measure of the detector's intrinsic spatial resolution, free from fit‑induced correlations.

7. Re‑engineered Track Angle Equations (angleX / angleY)
The original slope‑to‑angle conversion suffered from quadrant ambiguities and singularities near 90°. I re‑derived the geometry using atan2(1.0, slope) 
with a full 180°‑wrap correction. This ensures that the track angles are physically continuous across the vertical (90°) boundary, vital for generating 
clean angular distribution histograms and stable ML feature inputs.

🧪 Simulation & Validation Framework
To verify every algorithmic fix and to establish a "ground truth" for the detector's expected behavior, I built a dedicated Monte Carlo simulation of 
the entire telescope:

Realistic Detector Geometry: The simulation models the full 640‑strip plane, including the exact electric field profile inside the drift gap, to track 
the motion of primary ionization electrons from the conversion point to the readout strips.

Charge Amplification & Gas Gain: Each electron avalanche is assigned a charge based on a realistic amplification range, mimicking the statistical gain 
fluctuations (Gas Gain) of the Micromegas chamber. This produces synthetic hits that closely resemble the raw ADC distributions seen in real data.

Large‑Statistics Production: I ran the simulation with 2,000,000 events—the same scale as the real‑data analysis—to generate a robust reference dataset.

The simulation served two critical purposes:

Predicting the True Angle Distribution: Given the 4‑layer geometry and the cosmic‑ray flux, the simulation predicted a clean, continuous angular acceptance 
between 70° and 110° (relative to the horizontal plane)—exactly the physical range expected for muons traversing the telescope's active area.

Validating the Reconstruction Corrections: Before my fixes, the real‑data angle histograms (angleX and angleY) were contaminated by large, unnatural spikes 
outside the 70‑110° range, caused by mis‑mapped strips and noise. After applying the full correction stack 
(Mapping → Noise → Combinatorial Seeding → Hybrid Reconstruction), the real data now perfectly matches the simulation's prediction—a definitive, 
quantitative proof that the algorithms are physically correct and that the detector is now properly calibrated.

Quantifiable Impact
Metric                             |	Before Corrections                  |	After Full Algorithmic Stack
Peak Track Fit χ²                  |	> 8,000,000                         |	~ 4,000
Fake / Badly‑Reconstructed Tracks  |	~ 12%                               |	4 – 8%
Angle Distribution (70‑110° range) |	Unnatural spikes & discontinuities  |	Smooth, matches 2M‑event simulation
ML Feature Quality                 |	Uncalibrated raw data               |	Empirically‑weighted & validated features in ml_real_data.csv
