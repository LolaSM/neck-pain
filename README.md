# Neck Pain Classification Using Biomechanical Time-Series Data
This repository contains all the code, data, and processing scripts used in a study aiming to classify individuals with and without neck pain using biomechanical head motion data recorded during a VR-based functional task.

## Project Objective
The goal is to explore the potential of machine learning and deep learning models to classify subjects into two categories:
- **Group A**: pain-free (no neck pain)
- **Group B**: neck pain-affected

Data was captured using the HTC Vive Pro Eye system while participants performed the FarmDay task in virtual reality. From the 3D time-series data, temporal features such as velocity, acceleration, and displacement were extracted for analysis.

## Repository Structure

neck-pain/
 MATLAB/
   comparison_time_series_subjectA_subjectB.m
   patients_data_speed_acceleration.m
 RNN_LSTM.ipynb
 RNN_LSTM_dropout.ipynb
 comparison_performance_classification_models.ipynb
 tablaVectores_AB_speed_acceleration.csv
 tablagrupos_AB_raw.csv

## Folder and File Descriptions

### MATLAB/
Contains preprocessing scripts developed in **MATLAB R2022b** for:
- **`patients_data_speed_acceleration.m`**: Feature engineering -- calculation and extraction of statistical descriptors (mean, std, skw, kur) of velocity and acceleration from raw position data for each axis.
- **`comparison_time_series_subjectA_subjectB.m`**: time-series visualization for comparison between group A and B.

### Jupyter Notebooks 
Developed in Google Colab:
- **`comparison_performance_classification_models.ipynb`**:
  - Exploratory data analysis of the features with Histograms and Q-Q plots
  - Training and comparison of classical ML models (e.g., SVM, Random Forest, MLP) using manually extracted features of speed and acceleration and nested cross-validation.
- **`RNN_LSTM.ipynb`**: Baseline LSTM model for classifying neck pain based on raw time-series data. Requires GPU acceleration
- **`RNN_LSTM_dropout.ipynb`**: LSTM architecture adding a dropout layer. Requires GPU acceleration. Requires GPU acceleration

### CSV Files
- **`tablagrupos_AB_raw.csv`**: raw time-series data with class labels for each subject already preprocessed in MATLAB
- **`tablaVectores_AB_speed_acceleration.csv`**: statistical descriptors of velocity and acceleration vectors for each patient used for classical ML models and generated in MATLAB (`patients_data_speed_acceleration.m`)

## Dependencies
- Python 3.9.1
- TensorFlow 2.18.0, Keras 3.8.0
- Scikit-learn 1.6.1
- Seaborn 0.13.2
- NumPy 2.0.2, Pandas 2.2.2, Matplotlib 3.10.0
