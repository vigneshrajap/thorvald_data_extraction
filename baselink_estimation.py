from numpy import sin, cos, array
from math import atan2
from csv import reader, writer

# DCM matrix from b-frame to g-frame (NED)
def Cb_g(roll, pitch, yaw):
    return array([[cos(yaw)*cos(pitch), -sin(yaw)*cos(roll) + cos(yaw)*sin(pitch)*sin(roll), sin(yaw)*sin(roll) + cos(yaw)*sin(pitch)*cos(roll)],
                  [sin(yaw)*cos(pitch), cos(yaw)*cos(roll) + sin(yaw)*sin(pitch)*sin(roll), -cos(yaw)*sin(roll) + sin(yaw)*sin(pitch)*cos(roll)],
                  [-sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)]])

# Open files
with open('20191010_L4_N_slaloam/raw_gps.csv', 'r') as gnss_read_obj, \
     open('20191010_L4_N_slaloam/raw_imu.csv', 'r') as imu_read_obj, \
     open('20191010_L4_N_slaloam/baselink.csv', 'w') as base_write_obj:

    # Skip GNSS header
    gnss_data = reader(gnss_read_obj)
    gnss_line = next(gnss_data)

    # Skip IMU header
    imu_data = reader(imu_read_obj)
    imu_line = next(imu_data)


    # Estimating initial heading from GNSS using first 10sec of driving straight forward...
    # Time(Sec),Lat(deg),Long(deg),UTM_x(m),UTM_y(m)
    # 0.0000000000,59.6610704017,10.6748158033,6614856.7912911261,594364.6325763463
    # 10.6025481224,59.6610857330,10.6748542400,6614858.5529046170,594366.7548460397
    yaw_gnss = atan2((594366.755 - 594364.633), (6614858.553 - 6614856.791))

    # IMU heading
    # Time(Sec),Roll(rad),Pitch(rad),Yaw(rad),A_x(m/s2),A_y(m/s2),A_z(m/s2)
    # 0.000000,-0.050188,-0.106944,0.164520,1.036599,-0.506879,9.693512
    imu_first = next(imu_data)
    print(imu_first[3])
    yaw_imu = float(imu_first[3]) #0.164520
    magnetic_declination_radians = 0.0523599 
    
    # Estimating initial heading bias in IMU...
    yaw_bias = (yaw_gnss - yaw_imu) #+ magnetic_declination_radians

    # Eccentricities
    # GNSS relative to baselink in b-frame
    tx_gnss = 0.425         # meter
    ty_gnss = -0.620        # meter
    tz_gnss = 1.050         # meter
    T_gnss = array([[tx_gnss],
                    [ty_gnss],
                    [tz_gnss]])

    # Coordinate axis definitions
    NED2ENU = array([[0, 1, 0],
                     [1, 0, 0],
                     [0, 0, -1]])
    NWU2ENU = array([[0, -1, 0],
                     [1, 0, 0],
                     [0, 0, 1]])
    NWU2NED = array([[1, 0, 0],
                     [0, -1, 0],
                     [0, 0, -1]])

    # Synchronize observations
    for gnss_row in gnss_data:
        for imu_row in imu_data:
            if float(gnss_row[0]) <= float(imu_row[0]):
                sync_error = abs(float(gnss_row[0]) - float(imu_row[0]))    # sec
                N = float(gnss_row[3])                                      # meter
                E = float(gnss_row[4])                                      # meter
                U = 100                                                     # meter
                roll = float(imu_row[1])                                    # radian
                pitch = float(imu_row[2])                                   # radian
                yaw_imu = float(imu_row[3])                                 # radian

                P3_g = array([[E],
                              [N],
                              [U]])

                # IMU corrected heading (adding bias)
                yaw = yaw_imu + yaw_bias

                # Transformation from GNSS ENU g-frame to GNSS NWU b-frame
                P3_b = NWU2NED.T @ Cb_g(roll, pitch, yaw).T @ NED2ENU.T @ P3_g

                # Transformation from GNSS NWU b-frame to baselink NWU b-frame
                P1_b = P3_b - T_gnss

                # Transformation from baselink NWU b-frame to baselink ENU g-frame
                P1_g = NED2ENU @ Cb_g(roll, pitch, yaw) @ NWU2NED @ P1_b

                # Write to csv
                baselink = writer(base_write_obj)
                baselink.writerow([gnss_row[0], P1_g[1][0], P1_g[0][0], P1_g[2][0]])
                
                print(gnss_row[0], P1_g[1][0], P1_g[0][0], P1_g[2][0])
                
                break
