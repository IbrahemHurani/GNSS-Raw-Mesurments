
import csv
import math
from datetime import datetime, timezone, timedelta
import pandas as pd
import numpy as np
import simplekml
from cv2 import transform
from matplotlib import pyplot as plt
from navpy.core import navpy
from pyproj import Proj, transform
from scipy.optimize import least_squares
from simplekml import AltitudeMode, Kml

from gnssutils import EphemerisManager

WEEKSEC = 604800
LIGHTSPEED = 2.99792458e8
def forward_positioning_algorithm(sat_positions, pseudoranges, x0, b0):
    iterations = 0
    max_iterations = 10
    convergence_threshold = 1e-3

    while iterations < max_iterations:
        # Compute predicted pseudoranges
        r = np.linalg.norm(sat_positions - x0, axis=1)
        predicted_pseudoranges = r + b0

        # Ensure pseudoranges and predicted_pseudoranges match in length
        if len(pseudoranges) != len(predicted_pseudoranges):
            raise ValueError("Mismatch in array lengths: pseudoranges has length {}, predicted_pseudoranges has length {}".format(len(pseudoranges), len(predicted_pseudoranges)))

        # Compute residuals
        residuals = pseudoranges - predicted_pseudoranges

        # Check for convergence
        if np.max(np.abs(residuals)) < convergence_threshold:
            break

        # Formulate linearized system
        H = np.zeros((len(sat_positions), 4))
        for i in range(len(sat_positions)):
            H[i, :3] = -(sat_positions[i] - x0) / r[i]
            H[i, 3] = 1

        # Weighted least squares solution
        W = np.eye(len(sat_positions))
        x = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ residuals

        dx = x[:3]
        db = x[3]

        x0 = x0 + dx
        b0 = b0 + db

        iterations += 1

    return x0, b0

def calculate_satellite_position(ephemeris, transmit_time):
    mu = 3.986005e14
    OmegaDot_e = 7.2921151467e-5
    F = -4.442807633e-10
    sv_position = pd.DataFrame()
    sv_position['sv'] = ephemeris.index
    sv_position.set_index('sv', inplace=True)
    sv_position['t_k'] = transmit_time - ephemeris['t_oe']
    A = ephemeris['sqrtA'].pow(2)
    n_0 = np.sqrt(mu / A.pow(3))
    n = n_0 + ephemeris['deltaN']
    M_k = ephemeris['M_0'] + n * sv_position['t_k']
    E_k = M_k
    err = pd.Series(data=[1] * len(sv_position.index))
    i = 0
    while err.abs().min() > 1e-8 and i < 10:
        new_vals = M_k + ephemeris['e'] * np.sin(E_k)
        err = new_vals - E_k
        E_k = new_vals
        i += 1

    sinE_k = np.sin(E_k)
    cosE_k = np.cos(E_k)
    delT_r = F * ephemeris['e'].pow(ephemeris['sqrtA']) * sinE_k
    delT_oc = transmit_time - ephemeris['t_oc']
    sv_position['delT_sv'] = ephemeris['SVclockBias'] + ephemeris['SVclockDrift'] * delT_oc + ephemeris[
        'SVclockDriftRate'] * delT_oc.pow(2)

    v_k = np.arctan2(np.sqrt(1 - ephemeris['e'].pow(2)) * sinE_k, (cosE_k - ephemeris['e']))

    Phi_k = v_k + ephemeris['omega']

    sin2Phi_k = np.sin(2 * Phi_k)
    cos2Phi_k = np.cos(2 * Phi_k)

    du_k = ephemeris['C_us'] * sin2Phi_k + ephemeris['C_uc'] * cos2Phi_k
    dr_k = ephemeris['C_rs'] * sin2Phi_k + ephemeris['C_rc'] * cos2Phi_k
    di_k = ephemeris['C_is'] * sin2Phi_k + ephemeris['C_ic'] * cos2Phi_k

    u_k = Phi_k + du_k

    r_k = A * (1 - ephemeris['e'] * np.cos(E_k)) + dr_k

    i_k = ephemeris['i_0'] + di_k + ephemeris['IDOT'] * sv_position['t_k']

    x_k_prime = r_k * np.cos(u_k)
    y_k_prime = r_k * np.sin(u_k)

    Omega_k = ephemeris['Omega_0'] + (ephemeris['OmegaDot'] - OmegaDot_e) * sv_position['t_k'] - OmegaDot_e * \
              ephemeris[
                  't_oe']

    sv_position['x_k'] = x_k_prime * np.cos(Omega_k) - y_k_prime * np.cos(i_k) * np.sin(Omega_k)
    sv_position['y_k'] = x_k_prime * np.sin(Omega_k) + y_k_prime * np.cos(i_k) * np.cos(Omega_k)
    sv_position['z_k'] = y_k_prime * np.sin(i_k)
    return sv_position
def forward_positioning_algorithm(sat_positions, pseudoranges, x0, b0):
    iterations = 0
    max_iterations = 10
    convergence_threshold = 1e-3

    while iterations < max_iterations:
        # Compute predicted pseudoranges
        r = np.linalg.norm(sat_positions - x0, axis=1)
        predicted_pseudoranges = r + b0

        # Ensure pseudoranges and predicted_pseudoranges match in length
        if len(pseudoranges) != len(predicted_pseudoranges):
            raise ValueError("Mismatch in array lengths: pseudoranges has length {}, predicted_pseudoranges has length {}".format(len(pseudoranges), len(predicted_pseudoranges)))

        # Compute residuals
        residuals = pseudoranges - predicted_pseudoranges

        # Check for convergence
        if np.max(np.abs(residuals)) < convergence_threshold:
            break

        # Formulate linearized system
        H = np.zeros((len(sat_positions), 4))
        for i in range(len(sat_positions)):
            H[i, :3] = -(sat_positions[i] - x0) / r[i]
            H[i, 3] = 1

        # Weighted least squares solution
        W = np.eye(len(sat_positions))
        x = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ residuals

        dx = x[:3]
        db = x[3]

        x0 = x0 + dx
        b0 = b0 + db

        iterations += 1

    return x0, b0
def process(input_filepath):
    with open(input_filepath) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0][0] == '#':
                if 'Fix' in row[0]:
                    android_fixes = [row[1:]]
                elif 'Raw' in row[0]:
                    measurements = [row[1:]]
            else:
                if row[0] == 'Fix':
                    android_fixes.append(row[1:])
                elif row[0] == 'Raw':
                    measurements.append(row[1:])

    android_fixes = pd.DataFrame(android_fixes[1:], columns=android_fixes[0])
    measurements = pd.DataFrame(measurements[1:], columns=measurements[0])

    # Format satellite IDs
    measurements.loc[measurements['Svid'].str.len() == 1, 'Svid'] = '0' + measurements['Svid']
    measurements.loc[measurements['ConstellationType'] == '1', 'Constellation'] = 'G'
    measurements.loc[measurements['ConstellationType'] == '3', 'Constellation'] = 'R'
    measurements['SvName'] = measurements['Constellation'] + measurements['Svid']

    # Remove all non-GPS measurements
    measurements = measurements.loc[measurements['Constellation'] == 'G']

    # Convert columns to numeric representation
    measurements['Cn0DbHz'] = pd.to_numeric(measurements['Cn0DbHz'])
    measurements['TimeNanos'] = pd.to_numeric(measurements['TimeNanos'])
    measurements['FullBiasNanos'] = pd.to_numeric(measurements['FullBiasNanos'])
    measurements['ReceivedSvTimeNanos'] = pd.to_numeric(measurements['ReceivedSvTimeNanos'])
    measurements['PseudorangeRateMetersPerSecond'] = pd.to_numeric(measurements['PseudorangeRateMetersPerSecond'])
    measurements['ReceivedSvTimeUncertaintyNanos'] = pd.to_numeric(measurements['ReceivedSvTimeUncertaintyNanos'])

    # A few measurement values are not provided by all phones
    # We'll check for them and initialize them with zeros if missing
    if 'BiasNanos' in measurements.columns:
        measurements['BiasNanos'] = pd.to_numeric(measurements['BiasNanos'])
    else:
        measurements['BiasNanos'] = 0
    if 'TimeOffsetNanos' in measurements.columns:
        measurements['TimeOffsetNanos'] = pd.to_numeric(measurements['TimeOffsetNanos'])
    else:
        measurements['TimeOffsetNanos'] = 0

    # print(measurements.columns)

    measurements['GpsTimeNanos'] = measurements['TimeNanos'] - (
                measurements['FullBiasNanos'] - measurements['BiasNanos'])
    gpsepoch = datetime(1980, 1, 6, 0, 0, 0)
    measurements['UnixTime'] = pd.to_datetime(measurements['GpsTimeNanos'], utc=True, origin=gpsepoch)
    measurements['UnixTime'] = measurements['UnixTime']

    # Split data into measurement epochs
    measurements['Epoch'] = 0
    measurements.loc[
        measurements['UnixTime'] - measurements['UnixTime'].shift() > timedelta(milliseconds=200), 'Epoch'] = 1
    measurements['Epoch'] = measurements['Epoch'].cumsum()



    # This should account for rollovers since it uses a week number specific to each measurement

    measurements['tRxGnssNanos'] = measurements['TimeNanos'] + measurements['TimeOffsetNanos'] - (
                measurements['FullBiasNanos'].iloc[0] + measurements['BiasNanos'].iloc[0])
    measurements['GpsWeekNumber'] = np.floor(1e-9 * measurements['tRxGnssNanos'] / WEEKSEC)
    measurements['tRxSeconds'] = 1e-9 * measurements['tRxGnssNanos'] - WEEKSEC * measurements['GpsWeekNumber']
    measurements['tTxSeconds'] = 1e-9 * (measurements['ReceivedSvTimeNanos'] + measurements['TimeOffsetNanos'])
    # Calculate pseudorange in seconds
    measurements['prSeconds'] = measurements['tRxSeconds'] - measurements['tTxSeconds']

    # Conver to meters
    measurements['PrM'] = LIGHTSPEED * measurements['prSeconds']
    measurements['PrSigmaM'] = LIGHTSPEED * 1e-9 * measurements['ReceivedSvTimeUncertaintyNanos']

    manager = EphemerisManager()

    epoch = 1
    timestamp = measurements.iloc[epoch]['UnixTime'].to_pydatetime(warn=False)
    one_epoch = measurements.loc[(measurements['Epoch'] == epoch) & (measurements['prSeconds'] < 0.1)].drop_duplicates(
        subset='SvName')
    one_epoch.set_index('SvName', inplace=True)

    sats = one_epoch.index.unique().tolist()
    ephemeris = manager.get_ephemeris(timestamp, sats)
    b0 = 0
    x0 = np.array([0, 0, 0])

    List_forwadPos=[]
    ecef_list = []
    for epoch in measurements['Epoch'].unique():
        one_epoch = measurements.loc[(measurements['Epoch'] == epoch) & (measurements['prSeconds'] < 0.1)]
        one_epoch = one_epoch.drop_duplicates(subset='SvName').set_index('SvName')
        if len(one_epoch.index) > 4:
            timestamp = one_epoch.iloc[0]['UnixTime'].to_pydatetime(warn=False)
            sats = one_epoch.index.unique().tolist()
            ephemeris = manager.get_ephemeris(timestamp, sats)
            sv_position = calculate_satellite_position(ephemeris, one_epoch['tTxSeconds'])
            sv_position.index = sv_position.index.map(str)
            one_epoch = one_epoch.join(sv_position[['delT_sv']], how='left')
            one_epoch['PrM_corrected'] = one_epoch['PrM'] + LIGHTSPEED * one_epoch['delT_sv']
            xs = sv_position[['x_k', 'y_k', 'z_k']].to_numpy()
            pr = one_epoch['PrM'] + LIGHTSPEED * sv_position['delT_sv']
            pr = pr.to_numpy()
            x, b = forward_positioning_algorithm(xs, pr, x0, b0)
            List_forwadPos.append(x)


            for sv_Index in one_epoch.index:
                ecef_list.append({
                    "Time": timestamp.isoformat(),
                    "SatPRN (ID)": sv_Index,
                    "Pseudo-Range": one_epoch.at[sv_Index, 'PrM_corrected'],
                    "Cn0": one_epoch.at[sv_Index, 'Cn0DbHz'],
                    "Sat.x": sv_position.at[sv_Index, 'x_k'] if sv_Index in sv_position.index else np.nan,
                    "Sat.y": sv_position.at[sv_Index, 'y_k'] if sv_Index in sv_position.index else np.nan,
                    "Sat.z": sv_position.at[sv_Index, 'z_k'] if sv_Index in sv_position.index else np.nan

                })




    ecef_array = np.stack(List_forwadPos, axis=0)
    lla_array = np.stack(navpy.ecef2lla(ecef_array), axis=1)

    # Extract the first position as a reference for the NED transformation
    ref_lla = lla_array[0, :]
    ned_array = navpy.ecef2ned(ecef_array, ref_lla[0], ref_lla[1], ref_lla[2])

    # Convert back to Pandas and save to csv
    lla_df = pd.DataFrame(lla_array, columns=['Latitude', 'Longitude', 'Altitude'])
    ned_df = pd.DataFrame(ned_array, columns=['N', 'E', 'D'])
    lla_df.to_csv('calculated_postion.csv')
    android_fixes.to_csv('android_position.csv')

    # Plot
    plt.style.use('dark_background')
    plt.plot(ned_df['E'], ned_df['N'])
    plt.title('Position Offset from First Epoch')
    plt.xlabel("East (m)")
    plt.ylabel("North (m)")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

    return ecef_list

def csv_to_kml(csv_file, kml_file):
      kml = Kml()

      # Open the CSV file
      with open(csv_file, 'r') as file:
          reader = csv.DictReader(file)

          # Iterate over each row in the CSV
          for row in reader:
              # Get coordinates
              lat = float(row['LatitudeDegrees'])
              lon = float(row['LongitudeDegrees'])

              # Handle possible missing altitude data
              alt = row['AltitudeMeters']
              if alt:  # Check if altitude is not empty
                  alt = float(alt)
              else:
                  alt = 0.0  # Default altitude if missing

              # Create a point for this data
              pnt = kml.newpoint(name="Point", coords=[(lon, lat, alt)])
              pnt.altitudemode = AltitudeMode.absolute

      # Save the KML
      kml.save(kml_file)




if __name__ == '__main__':
  ecef_List=process('C:\\Users\\user\Desktop\gnss-analysis-main\data\sample\walking.txt')
  csv_FromList = pd.DataFrame(ecef_List)
  csv_FromList.to_csv("Measurements_Output.csv", index=False)
  #####Q3####
  file_measurements_output='android_position.csv'
  csv_to_kml(file_measurements_output,'output_kml.kml')












