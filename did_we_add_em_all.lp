\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  7.920951 gen[1] + 23.269494 gen[2]
Subject To
 F_eq[1,2]: - corr_flow[1,2] - 18.93939393939394 bus_angle[1]
   + 18.93939393939394 bus_angle[2] = 0
 F_eq[1,5]: - corr_flow[1,5] - 20.32520325203252 bus_angle[1]
   + 20.32520325203252 bus_angle[5] = 0
 F_eq[2,3]: - corr_flow[2,3] - 22.8310502283105 bus_angle[2]
   + 22.8310502283105 bus_angle[3] = 0
 F_eq[2,4]: - corr_flow[2,4] - 29.41176470588235 bus_angle[2]
   + 29.41176470588235 bus_angle[4] = 0
 F_eq[2,5]: - corr_flow[2,5] - 28.90173410404624 bus_angle[2]
   + 28.90173410404624 bus_angle[5] = 0
 F_eq[3,4]: - corr_flow[3,4] - 78.125 bus_angle[3] + 78.125 bus_angle[4]
   = 0
 F_eq[4,5]: - corr_flow[4,5] - 23.74732842555212 bus_angle[4]
   + 23.74732842555212 bus_angle[5] = 0
 F_eq[4,7]: - corr_flow[4,7] - 4.781943381790359 bus_angle[4]
   + 4.781943381790359 bus_angle[7] = 0
 F_eq[4,9]: - corr_flow[4,9] - 1.797979071523607 bus_angle[4]
   + 1.797979071523607 bus_angle[9] = 0
 F_eq[5,6]: - corr_flow[5,6] - 3.967939052456154 bus_angle[5]
   + 3.967939052456154 bus_angle[6] = 0
 F_eq[6,11]: - corr_flow[6,11] - 5.027652086475616 bus_angle[6]
   + 5.027652086475616 bus_angle[11] = 0
 F_eq[6,12]: - corr_flow[6,12] - 3.909151323247723 bus_angle[6]
   + 3.909151323247723 bus_angle[12] = 0
 F_eq[6,13]: - corr_flow[6,13] - 7.676364473785216 bus_angle[6]
   + 7.676364473785216 bus_angle[13] = 0
 F_eq[7,8]: - corr_flow[7,8] - 5.676979846721544 bus_angle[7]
   + 5.676979846721544 bus_angle[8] = 0
 F_eq[7,9]: - corr_flow[7,9] - 9.090082719752751 bus_angle[7]
   + 9.090082719752751 bus_angle[9] = 0
 F_eq[9,10]: - corr_flow[9,10] - 11.83431952662722 bus_angle[9]
   + 11.83431952662722 bus_angle[10] = 0
 F_eq[9,14]: - corr_flow[9,14] - 3.698498409645684 bus_angle[9]
   + 3.698498409645684 bus_angle[14] = 0
 F_eq[10,11]: - corr_flow[10,11] - 5.206435153850159 bus_angle[10]
   + 5.206435153850159 bus_angle[11] = 0
 F_eq[12,13]: - corr_flow[12,13] - 5.003001801080648 bus_angle[12]
   + 5.003001801080648 bus_angle[13] = 0
 F_eq[13,14]: - corr_flow[13,14] - 2.873398080570082 bus_angle[13]
   + 2.873398080570082 bus_angle[14] = 0
 R20: - corr_flow[1,2] - corr_flow[1,5] + gen[1] = 0
 R21: corr_flow[1,2] - corr_flow[2,3] - corr_flow[2,4] - corr_flow[2,5]
   + gen[2] = 8.68
 R22: corr_flow[2,3] - corr_flow[3,4] + gen[3] = 37.68
 R23: corr_flow[2,4] + corr_flow[3,4] - corr_flow[4,5] - corr_flow[4,7]
   - corr_flow[4,9] + gen[4] = 19.12
 R24: corr_flow[1,5] + corr_flow[2,5] + corr_flow[4,5] - corr_flow[5,6]
   + gen[5] = 3.04
 R25: corr_flow[5,6] - corr_flow[6,11] - corr_flow[6,12] - corr_flow[6,13]
   + gen[6] = 4.48
 R26: corr_flow[4,7] - corr_flow[7,8] - corr_flow[7,9] + gen[7] = 0
 R27: corr_flow[7,8] + gen[8] = 0
 R28: corr_flow[4,9] + corr_flow[7,9] - corr_flow[9,10] - corr_flow[9,14]
   + gen[9] = 11.8
 R29: corr_flow[9,10] - corr_flow[10,11] + gen[10] = 3.6
 R30: corr_flow[6,11] + corr_flow[10,11] + gen[11] = 1.4
 R31: corr_flow[6,12] - corr_flow[12,13] + gen[12] = 2.44
 R32: corr_flow[6,13] + corr_flow[12,13] - corr_flow[13,14] + gen[13]
   = 5.4
 R33: corr_flow[9,14] + corr_flow[13,14] + gen[14] = 5.96
 R34: - bus_angle[9] + bus_angle[11] <= 93.461557
 R35: bus_angle[9] - bus_angle[11] <= 93.461557
 R36: bus_angle[11] - bus_angle[14] <= 110.294555
 R37: - bus_angle[11] + bus_angle[14] <= 110.294555
 R38: - bus_angle[1] + bus_angle[2] <= 39.534152
 R39: bus_angle[1] - bus_angle[2] <= 39.534152
 R40: bus_angle[4] - bus_angle[5] <= 112.217364
 R41: - bus_angle[4] + bus_angle[5] <= 112.217364
 R42: bus_angle[2] - bus_angle[4] <= 38.656352
 R43: - bus_angle[2] + bus_angle[4] <= 38.656352
 R44: bus_angle[2] - bus_angle[5] <= 39.169004
 R45: - bus_angle[2] + bus_angle[5] <= 39.169004
 R46: bus_angle[4] - bus_angle[6] <= 59.85391299999999
 R47: - bus_angle[4] + bus_angle[6] <= 59.85391299999999
 R48: bus_angle[4] - bus_angle[13] <= 114.742026
 R49: - bus_angle[4] + bus_angle[13] <= 114.742026
 R50: - bus_angle[6] + bus_angle[11] <= 113.880305
 R51: bus_angle[6] - bus_angle[11] <= 113.880305
 R52: bus_angle[11] - bus_angle[13] <= 59.75423300000002
 R53: - bus_angle[11] + bus_angle[13] <= 59.75423300000002
 R54: bus_angle[4] - bus_angle[9] <= 41.19514699999999
 R55: - bus_angle[4] + bus_angle[9] <= 41.19514699999999
 R56: - bus_angle[6] + bus_angle[12] <= 36.95195699999999
 R57: bus_angle[6] - bus_angle[12] <= 36.95195699999999
 R58: bus_angle[4] - bus_angle[14] <= 95.79415299999999
 R59: - bus_angle[4] + bus_angle[14] <= 95.79415299999999
Bounds
 corr_flow[1,2] free
 corr_flow[1,5] free
 corr_flow[2,3] free
 corr_flow[2,4] free
 corr_flow[2,5] free
 corr_flow[3,4] free
 corr_flow[4,5] free
 corr_flow[4,7] free
 corr_flow[4,9] free
 corr_flow[5,6] free
 corr_flow[6,11] free
 corr_flow[6,12] free
 corr_flow[6,13] free
 corr_flow[7,8] free
 corr_flow[7,9] free
 corr_flow[9,10] free
 corr_flow[9,14] free
 corr_flow[10,11] free
 corr_flow[12,13] free
 corr_flow[13,14] free
 gen[1] <= 340
 gen[2] <= 59
 gen[3] = 0
 gen[4] = 0
 gen[5] = 0
 gen[6] = 0
 gen[7] = 0
 gen[8] = 0
 gen[9] = 0
 gen[10] = 0
 gen[11] = 0
 gen[12] = 0
 gen[13] = 0
 gen[14] = 0
 -30 <= bus_angle[1] <= 30
 -30 <= bus_angle[2] <= 30
 -30 <= bus_angle[3] <= 30
 -30 <= bus_angle[4] <= 30
 -30 <= bus_angle[5] <= 30
 -30 <= bus_angle[6] <= 30
 -30 <= bus_angle[7] <= 30
 -30 <= bus_angle[8] <= 30
 -30 <= bus_angle[9] <= 30
 -30 <= bus_angle[10] <= 30
 -30 <= bus_angle[11] <= 30
 -30 <= bus_angle[12] <= 30
 -30 <= bus_angle[13] <= 30
 -30 <= bus_angle[14] <= 30
End
