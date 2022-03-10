# 3D Sonic Post Processing

Scripts to process high-frequency 3D sonic wind and temperature measurements into wind statistics, turbulent fluxes and derived quantities.

## Application

Run the script `process_sonic_55m_110m_sample` to post-process the Gill sonics (at 55m and 110m), which calls the other functions. You need to have the OneDAS connector, otherwise change the _filename_ path.

Run the script `process_sonic_25m_sample` to post-process the 25-m Thies sonic, which calls the other functions. You need to have the OneDAS connector, otherwise change the _filename_ path.

## Hourly processing to OneDAS (from Remote PC)

The scripts `process_sonic_55m_110m_hourly` and `process_sonic_25m_hourly` have the correct paths to load data from Z and save hourly post-processed ASCII files in the folder loaded by OneDAS.
