# Acoustic Gunshot Detection & TDOA Localization

This project implements an acoustic gunshot detection system using Machine Learning (CNN & Bagged Trees) and simulates source localization using Time Difference of Arrival (TDOA) triangulation in Simulink.

The system is designed to classify audio as **Gunshot** or **Non-Gunshot** and, upon detection, estimate the source coordinates using a simulated three-node sensor array.

---

##  Project Overview

### 1. Machine Learning Training

Two different machine learning approaches are used to distinguish **Gunshot** sounds from **Non-Gunshot** background noise.  
Both training scripts include **custom data-balancing techniques** to compensate for the short duration of gunshot events.

---

#### 🔹 Ensemble Classifier (Bagged Trees)

- **Script:** `xgboost_train.m`
- **Features:** MFCCs (Mel-Frequency Cepstral Coefficients)
- **Window Length:** 2 seconds
- **Architecture:** Bagged Ensemble of 100 Decision Trees
- **Output Model:** `gunshotModel.mat`
- **Usage:** Loaded into the Simulink *Predict* block

---

#### 🔹 Deep Learning Classifier (CNN)

- **Script:** `CNN_train.m`
- **Features:** Log-Mel Spectrograms
- **Architecture:**  
  - 2 Convolution Blocks  
  - Batch Normalization  
  - Max Pooling
- **Output Model:** `gunshotCNNModel.mat`
- **Usage:** Loaded via workspace variable in Simulink

---

### 2. Simulink Simulation Models

The repository includes **three Simulink models**, each serving a distinct testing purpose:

| Model Name | Input Source | Description |
|-----------|-------------|-------------|
| **SingleMCU_Test** | Multimedia File | **ML Validation Only.** Directly feeds an audio file into the ML classifier to verify detection accuracy without TDOA. |
| **SingleFileTest_TDOA** | Multimedia File | **TDOA Simulation.** Applies artificial delays to simulate reception at three MCU nodes and estimates source coordinates. |
| **RealTimeMic_TDOA** | Laptop Microphone | **Live Simulation.** Captures real-time audio, applies simulated propagation delays, and performs real-time triangulation. |

---

## CRITICAL: Run This First

**You must run the training scripts before opening Simulink.**

The Simulink models depend on specific workspace variables and `.mat` files to function. If you do not run the training scripts first, the **Predict Block** and **Classification** subsystems will fail to load.

1. **For Ensemble / MFCC Model**  
   Run `xgboost_train.m`  
   → Generates `gunshotModel.mat`
   → `trainedmodel` variable used in Predict Block

2. **For CNN / Spectrogram Model**  
   Run `CNN_train.m`  
   → Generates `gunshotCNNModel.mat`
   → which used in CNN-Block
   
4. **Data Agumentation**
   Run `data_agument.mat`
   → Generates augmented data on sound by adding noise to file also it is incorporated with doubling data Dateset U can comment it if not needed

5. **Delay**
    > Since in this simulation we have to consider different delay for three MCU node which are given manualy in simulink model so carefully select the delay which is main component of TDOA also same for MIC Positions
---

## 🚀 Usage Instructions

### Prerequisites

- MATLAB R2023a or later (Recommended)
- Audio Toolbox
- Deep Learning Toolbox
- DSP System Toolbox
- Simulink

---

### Dataset Structure

Ensure the files all on same folder
