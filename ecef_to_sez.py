# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
# Converts ECEF coordinate vector to SEZ coordinate vector
# 
# Parameters:
#  o_x_km: X distance value of origin in ECEF frame in kilometers
#  o_y_km: Y distance value of origin in ECEF frame in kilometers
#  o_z_km: Z distance value of origin in ECEF frame in kilometers
#  x_km: X distance value in ECEF frame in kilometers
#  y_km: Y distance value in ECEF frame in kilometers
#  z_km: Z distance value in ECEF frame in kilometers
#  
# Output:
#  Print the SEZ coordinates
#
# Written by Matthew Moore
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.

# "constants"
R_E_KM = 6378.1363
E_E = 0.081819221456

# import Python modules
# e.g., import math # math module
import math # math module
import sys # argv

# helper functions

## matrix_times_vector
##
def matrix_times_vector(mat, vec, out_vec):
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            out_vec[i] += mat[i][j] * vec[j]
    return out_vec

## calc_denom
##
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0-ecc**2.0 * math.sin(lat_rad)**2.0)
            

# initialize script arguments
rECEF_sez = [0, 0, 0]
Ry = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
Rz = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

# parse script arguments
if len(sys.argv)==7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km = float(sys.argv[3])
    x_km = float(sys.argv[4])
    y_km = float(sys.argv[5])
    z_km = float(sys.argv[6])
else:
    print('Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km ')
    exit()

# write script below this line

# determine llh from origin ECEF vector
lon_rad = math.atan2(o_y_km,o_x_km)
lon_deg = lon_rad*180.0/math.pi

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1
lat_deg = lat_rad * 180.0 / math.pi  

# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

# initialize and remove SEZ origin from ECEF vector
rECEF = [x_km - o_x_km, y_km - o_y_km, z_km - o_z_km]

# initialize calculation matrices
Ry = [[math.sin(lat_rad), 0, -math.cos(lat_rad)],
      [0, 1, 0],
      [math.cos(lat_rad), 0, math.sin(lat_rad)]]

Rz = [[math.cos(lon_rad), math.sin(lon_rad), 0],
      [-math.sin(lon_rad), math.cos(lon_rad), 0],
      [0, 0, 1]]

# use inverse rotation matrices to calculate SEZ vector
RzrECEF = [0, 0, 0]
rSEZ = [0, 0, 0]
RzrECEF = matrix_times_vector(Rz, rECEF, RzrECEF)
rSEZ = matrix_times_vector(Ry, RzrECEF, rSEZ)

# place SEZ vector into individual variables
s_km = rSEZ[0]
e_km = rSEZ[1]
z_km = rSEZ[2]

# output variables individually with no clarifying user interface
print(s_km)
print(e_km)
print(z_km)
