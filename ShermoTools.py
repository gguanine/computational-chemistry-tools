import cclib
import os
import subprocess
import platform
import re
import csv

# Maybe only support DFT in this version
# Shermo developed by Tian Lu et. al. (10.1016/j.comptc.2021.113249)

DIR_OPT = "opt"
DIR_FREQ = "freq"
DIR_ENERGY = "energy"

# Define path to Shermo here
system_platfrom = platform.system()
if  system_platfrom == "Windows":
    shermo_path = "Shermo.exe"
elif system_platfrom == "Linux":
    shermo_path = "Shermo"
else:
    shermo_path = ""


ThermoScaler = {
# Here use part of the data collected by @liyuanhe21 (https://github.com/liyuanhe211/Collection_of_Frequency_Scale_Factors)
# Truhlar: http://comp.chem.umn.edu/freqscale/version3b2.htm
# Keys should be in lowercase
# Key1|Method    Key2|BasisSet    Value[0]|ZPE    Value[1]|H(T)-H(0)    Value[2]|S(T)    Value[4]|ref
 'b1b95': {'6-31+g(d,p)': [0.976, 0.9838, 0.9874, '10.1021/jp073974n'],
           '6-31g(2df,p)': [0.9762, 0.9714, 0.9767, '10.1021/jp073974n'],
           '6-31g(d)': [0.9716, 0.9788, 0.9816, '10.1021/jp073974n']},
 'b1lyp': {'6-31+g(d,p)': [0.9805, 0.9977, 1.0011, '10.1021/jp073974n'],
           '6-311+g(2df,p)': [0.984, 0.986, 0.9906, '10.1021/jp073974n'],
           '6-311+g(d,p)': [0.9838, 1.0017, 1.0074, '10.1021/jp073974n'],
           '6-31g(2df,p)': [0.9802, 0.983, 0.9869, '10.1021/jp073974n'],
           '6-31g(d)': [0.976, 0.992, 0.9947, '10.1021/jp073974n']},
 'b3lyp': {'6-31+g(d)': [0.9829, 1.0053, 1.008, '10.1021/jp073974n'],
           '6-31+g(d,p)': [0.9857, 1.0062, 1.0099, '10.1021/jp073974n'],
           '6-311++g(3df,3pd)': [0.9876, 0.9904, 0.995, '10.1021/jp073974n'],
           '6-311+g(2d,p)': [0.9898, 1.004, 1.0094, '10.1021/jp073974n'],
           '6-311+g(2df,p)': [0.9889, 0.9938, 0.9984, '10.1021/jp073974n'],
           '6-311+g(3d,p)': [0.9891, 0.9991, 1.0041, '10.1021/jp073974n'],
           '6-311+g(3df,2p)': [0.9877, 0.9915, 0.9963, '10.1021/jp073974n'],
           '6-311+g(3df,3pd)': [0.9876, 0.9904, 0.995, '10.1021/jp073974n'],
           '6-311+g(3df,p)': [0.9883, 0.9919, 0.997, '10.1021/jp073974n'],
           '6-311+g(d)': [0.9882, 1.0124, 1.0174, '10.1021/jp073974n'],
           '6-311+g(d,p)': [0.9887, 1.0102, 1.0161, '10.1021/jp073974n'],
           '6-311g(d)': [0.9882, 1.009, 1.0123, '10.1021/jp073974n'],
           '6-311g(d,p)': [0.9888, 1.0062, 1.0104, '10.1021/jp073974n'],
           '6-31g(2df,p)': [0.9853, 0.9909, 0.9946, '10.1021/jp073974n'],
           '6-31g(d)': [0.9813, 1.0004, 1.0029, '10.1021/jp073974n'],
           '6-31g(d,p)': [0.9838, 1.0003, 1.0033, '10.1021/jp073974n'],
           'def2-tzvp': [0.9850, 1., 1., 'truhlar'],
           'aug-cc-pvdz': [0.9787, 0.9789, 0.9486, '10.1021/jp048233q'],
           'aug-cc-pvqz': [0.9884, 0.9859, 0.9673, '10.1021/jp048233q'],
           'aug-cc-pvtz': [0.9867, 0.9867, 0.9538, '10.1021/jp048233q'],
           'cc-pvdz': [0.9689, 0.9784, 0.9475, '10.1021/jp048233q'],
           'cc-pvqz': [0.9854, 0.9872, 0.9586, '10.1021/jp048233q'],
           'cc-pvtz': [0.9764, 0.9854, 0.9576, '10.1021/jp048233q'],
           'cc-pvtz+d': [0.9886, 0.9926, 0.997, '10.1021/acs.jctc.6b00554'],
           'pc-0': [0.987, 0.9687, 0.9438, '10.1002/jcc.23073'],
           'pc-1': [0.988, 1.014, 1.015, '10.1002/jcc.23073'],
           'pc-2': [0.9869, 1.017, 1.015, '10.1002/jcc.23073'],
           'pc-3': [0.9876, 1.017, 1.015, '10.1002/jcc.23073'],
           'pc-4': [0.9877, 1.017, 1.015, '10.1002/jcc.23073']},
 'b3p86': {'6-31+g(d,p)': [0.9809, 0.9927, 0.9974, '10.1021/jp073974n'],
           '6-311+g(2df,p)': [0.9852, 0.9811, 0.9867, '10.1021/jp073974n'],
           '6-311+g(d,p)': [0.9845, 0.9963, 1.003, '10.1021/jp073974n'],
           '6-31g(2df,p)': [0.9814, 0.9795, 0.9846, '10.1021/jp073974n'],
           '6-31g(d)': [0.9768, 0.9881, 0.9919, '10.1021/jp073974n'],
           'pc-0': [0.9801, 0.9794, 0.9558, '10.1002/jcc.23073'],
           'pc-1': [0.9819, 1.03, 1.032, '10.1002/jcc.23073'],
           'pc-2': [0.9834, 1.03, 1.028, '10.1002/jcc.23073'],
           'pc-3': [0.9844, 1.03, 1.028, '10.1002/jcc.23073'],
           'pc-4': [0.9845, 1.03, 1.029, '10.1002/jcc.23073']},
 'b3pw91': {'6-31+g(d,p)': [0.9819, 0.9947, 0.9994, '10.1021/jp073974n'],
            '6-311+g(2df,p)': [0.9865, 0.983, 0.9886, '10.1021/jp073974n'],
            '6-311+g(d,p)': [0.9858, 0.9978, 1.004, '10.1021/jp073974n'],
            '6-31g(2df,p)': [0.9825, 0.9815, 0.9866, '10.1021/jp073974n'],
            '6-31g(d)': [0.978, 0.9899, 0.9937, '10.1021/jp073974n'],
            'pc-0': [0.9835, 0.9766, 0.9531, '10.1002/jcc.23073'],
            'pc-1': [0.9832, 1.028, 1.029, '10.1002/jcc.23073'],
            'pc-2': [0.9851, 1.026, 1.025, '10.1002/jcc.23073'],
            'pc-3': [0.9861, 1.026, 1.025, '10.1002/jcc.23073'],
            'pc-4': [0.9861, 1.026, 1.025, '10.1002/jcc.23073']},
 'b97d': {'pc-0': [1.023, 0.9114, 0.8848, '10.1002/jcc.23073'],
           'pc-1': [1.012, 0.949, 0.9463, '10.1002/jcc.23073'],
           'pc-2': [1.011, 0.9528, 0.9481, '10.1002/jcc.23073'],
           'pc-3': [1.011, 0.9536, 0.9492, '10.1002/jcc.23073'],
           'pc-4': [1.012, 0.9427, 0.909, '10.1002/jcc.23073']},
 'b971': {'6-31+g(d,p)': [0.9859, 1.0037, 1.0074, '10.1021/jp073974n'],
          '6-311+g(2df,p)': [0.9899, 0.9922, 0.9968, '10.1021/jp073974n'],
          '6-311+g(d,p)': [0.9893, 1.0072, 1.0134, '10.1021/jp073974n'],
          '6-31g(2df,p)': [0.9859, 0.9899, 0.9939, '10.1021/jp073974n'],
          '6-31g(d)': [0.9817, 0.9989, 1.0017, '10.1021/jp073974n']},
 'b972': {'6-31+g(d,p)': [0.976, 0.9895, 0.9943, '10.1021/jp073974n'],
          '6-311+g(2df,p)': [0.9809, 0.9773, 0.9827, '10.1021/jp073974n'],
          '6-311+g(d,p)': [0.9799, 0.9922, 0.9983, '10.1021/jp073974n'],
          '6-31g(2df,p)': [0.9768, 0.9763, 0.981, '10.1021/jp073974n'],
          '6-31g(d)': [0.9719, 0.9847, 0.9883, '10.1021/jp073974n']},
 'b98': {'6-31+g(d,p)': [0.985, 1.0016, 1.0055, '10.1021/jp073974n'],
         '6-311+g(2df,p)': [0.9886, 0.9899, 0.9947, '10.1021/jp073974n'],
         '6-311+g(d,p)': [0.9884, 1.0046, 1.0106, '10.1021/jp073974n'],
         '6-31g(2df,p)': [0.9849, 0.9879, 0.9921, '10.1021/jp073974n'],
         '6-31g(d)': [0.9809, 0.9967, 0.9997, '10.1021/jp073974n']},
 'bhandh': {'6-31+g(d,p)': [0.9562, 0.9335, 0.9391, '10.1021/jp073974n'],
            '6-311+g(2df,p)': [0.962, 0.9247, 0.931, '10.1021/jp073974n'],
            '6-311+g(d,p)': [0.9607, 0.9352, 0.9412, '10.1021/jp073974n'],
            '6-31g(2df,p)': [0.9568, 0.9234, 0.9291, '10.1021/jp073974n'],
            '6-31g(d)': [0.95, 0.9276, 0.9322, '10.1021/jp073974n']},
 'bhandhlyp': {'6-31+g(d,p)': [0.9498, 0.9453, 0.951, '10.1021/jp073974n'],
               '6-311+g(2df,p)': [0.9547, 0.9361, 0.9425, '10.1021/jp073974n'],
               '6-311+g(d,p)': [0.954, 0.9478, 0.9541, '10.1021/jp073974n'],
               '6-31g(2df,p)': [0.9506, 0.9341, 0.9399, '10.1021/jp073974n'],
               '6-31g(d)': [0.9446, 0.94, 0.9449, '10.1021/jp073974n']},
 'blyp': {'6-31+g(d,p)': [1.0169, 1.0703, 1.0749, '10.1021/jp073974n'],
          '6-311+g(2df,p)': [1.0186, 1.0538, 1.0571, '10.1021/jp073974n'],
          '6-311+g(d,p)': [1.0189, 1.0766, 1.087, '10.1021/jp073974n'],
          '6-31g(2df,p)': [1.0158, 1.0504, 1.0527, '10.1021/jp073974n'],
          '6-31g(d)': [1.0135, 1.0648, 1.0683, '10.1021/jp073974n'],
          'pc-0': [1.018, 0.9094, 0.8837, '10.1002/jcc.23073'],
          'pc-1': [1.02, 0.9364, 0.932, '10.1002/jcc.23073'],
          'pc-2': [1.016, 0.9394, 0.9331, '10.1002/jcc.23073'],
          'pc-3': [1.017, 0.9401, 0.9337, '10.1002/jcc.23073'],
          'pc-4': [1.017, 0.9404, 0.9341, '10.1002/jcc.23073']},
 'bmk': {'6-31+g(d)': [0.9733, 0.9721, 0.9772, '10.1021/jp073974n'],
         '6-31+g(d,p)': [0.9773, 0.9728, 0.9781, '10.1021/jp073974n'],
         '6-311++g(3df,3pd)': [0.9779, 0.9613, 0.9679, '10.1021/jp073974n'],
         '6-311+g(2d,p)': [0.9788, 0.971, 0.9778, '10.1021/jp073974n'],
         '6-311+g(2df,p)': [0.9787, 0.9644, 0.9709, '10.1021/jp073974n'],
         '6-311+g(3d,p)': [0.9776, 0.9667, 0.9729, '10.1021/jp073974n'],
         '6-311+g(3df,2p)': [0.9772, 0.9619, 0.9684, '10.1021/jp073974n'],
         '6-311+g(3df,3pd)': [0.9779, 0.9614, 0.9679, '10.1021/jp073974n'],
         '6-311+g(3df,p)': [0.9784, 0.9623, 0.9689, '10.1021/jp073974n'],
         '6-311+g(d)': [0.9773, 0.975, 0.9817, '10.1021/jp073974n'],
         '6-311+g(d,p)': [0.9794, 0.974, 0.9813, '10.1021/jp073974n'],
         '6-311g(d)': [0.9767, 0.9727, 0.9787, '10.1021/jp073974n'],
         '6-311g(d,p)': [0.9788, 0.9713, 0.9781, '10.1021/jp073974n'],
         '6-31g(2df,p)': [0.9752, 0.961, 0.968, '10.1021/jp073974n'],
         '6-31g(d)': [0.9709, 0.9679, 0.9731, '10.1021/jp073974n'],
         '6-31g(d,p)': [0.975, 0.9682, 0.9739, '10.1021/jp073974n'],
         'tzv2p': [0.976, 0.9722, 0.9798, '10.1021/jp073974n'],
         'aug-cc-pvdz': [0.9838, 0.9988, 1.0121, '10.1021/jp073974n'],
         'aug-cc-pvqz': [0.9777, 0.9657, 0.9741, '10.1021/jp073974n'],
         'aug-cc-pvtz': [0.9802, 0.9698, 0.9778, '10.1021/jp073974n'],
         'cc-pvdz': [0.9813, 0.9771, 0.9826, '10.1021/jp073974n'],
         'cc-pvqz': [0.977, 0.9643, 0.9724, '10.1021/jp073974n'],
         'cc-pvtz': [0.9792, 0.9675, 0.9758, '10.1021/jp073974n']},
 'bp86': {'6-31+g(d,p)': [1.0155, 1.0546, 1.0602, '10.1021/jp073974n'],
          '6-311+g(2df,p)': [1.0185, 1.0387, 1.043, '10.1021/jp073974n'],
          '6-311+g(d,p)': [1.0183, 1.0603, 1.0726, '10.1021/jp073974n'],
          '6-31g(2df,p)': [1.015, 1.0368, 1.0409, '10.1021/jp073974n'],
          '6-31g(d)': [1.0121, 1.0502, 1.0544, '10.1021/jp073974n'],
          'pc-0': [1.016, 0.9186, 0.8949, '10.1002/jcc.23073'],
          'pc-1': [1.017, 0.9493, 0.9468, '10.1002/jcc.23073'],
          'pc-2': [1.017, 0.9508, 0.9458, '10.1002/jcc.23073'],
          'pc-3': [1.017, 0.9509, 0.9458, '10.1002/jcc.23073'],
          'pc-4': [1.017, 0.9511, 0.9462, '10.1002/jcc.23073']},
 'ccsd': {'6-31+g(d,p)': [0.9686, 1.0043, 1.0145, '10.1021/jp073974n'],
          '6-311+g(d,p)': [0.9795, 0.9918, 0.9998, '10.1021/jp073974n'],
          '6-31g(d)': [0.9758, 0.9936, 1.0009, '10.1021/jp073974n']},
 'ccsd(t)': {'6-31+g(d,p)': [0.9779, 1.0391, 1.0555, '10.1021/jp073974n'],
             '6-311+g(d,p)': [0.9897, 1.0244, 1.0351, '10.1021/jp073974n'],
             '6-31g(d)': [0.9851, 1.0247, 1.0376, '10.1021/jp073974n']},
 'hcth147': {'6-31+g(d,p)': [0.9999, 1.0372, 1.042, '10.1021/jp073974n'],
             '6-311+g(2df,p)': [1.0037, 1.0199, 1.0233, '10.1021/jp073974n'],
             '6-311+g(d,p)': [1.0028, 1.0404, 1.0486, '10.1021/jp073974n'],
             '6-31g(2df,p)': [0.9998, 1.0193, 1.022, '10.1021/jp073974n'],
             '6-31g(d)': [0.9964, 1.0325, 1.0357, '10.1021/jp073974n']},
 'hcth407': {'6-31+g(d,p)': [0.995, 1.0292, 1.0342, '10.1021/jp073974n'],
             '6-311+g(2df,p)': [0.9993, 1.0124, 1.0168, '10.1021/jp073974n'],
             '6-311+g(d,p)': [0.9981, 1.0322, 1.0399, '10.1021/jp073974n'],
             '6-31g(2df,p)': [0.9951, 1.0118, 1.0149, '10.1021/jp073974n'],
             '6-31g(d)': [0.9911, 1.024, 1.0273, '10.1021/jp073974n']},
 'hcth93': {'6-31+g(d,p)': [0.9991, 1.0357, 1.041, '10.1021/jp073974n'],
            '6-311+g(2df,p)': [1.0033, 1.0185, 1.0224, '10.1021/jp073974n'],
            '6-311+g(d,p)': [1.0022, 1.0383, 1.0466, '10.1021/jp073974n'],
            '6-31g(2df,p)': [0.9992, 1.0182, 1.0214, '10.1021/jp073974n'],
            '6-31g(d)': [0.9957, 1.0311, 1.0347, '10.1021/jp073974n'],
            'pc-0': [1.002, 0.933, 0.9108, '10.1002/jcc.23073'],
            'pc-1': [1.0, 0.9806, 0.9811, '10.1002/jcc.23073'],
            'pc-2': [1.002, 0.9788, 0.9765, '10.1002/jcc.23073'],
            'pc-3': [1.003, 0.9789, 0.9765, '10.1002/jcc.23073'],
            'pc-4': [1.003, 0.9787, 0.9763, '10.1002/jcc.23073']},
 'hf': {'3-21g': [0.9207, 0.9444, 0.9666, '10.1021/jp960976r'],
        '6-31+g(d)': [0.9154, 0.8943, 0.9024, '10.1021/jp073974n'],
        '6-31+g(d,p)': [0.92, 0.8955, 0.9035, '10.1021/jp073974n'],
        '6-311++g(3df,3pd)': [0.9265, 0.8864, 0.8938, '10.1021/jp073974n'],
        '6-311+g(2d,p)': [0.9273, 0.8964, 0.9044, '10.1021/jp073974n'],
        '6-311+g(2df,p)': [0.9268, 0.8885, 0.8959, '10.1021/jp073974n'],
        '6-311+g(3d,p)': [0.9268, 0.8938, 0.9013, '10.1021/jp073974n'],
        '6-311+g(3df,2p)': [0.926, 0.8875, 0.8949, '10.1021/jp073974n'],
        '6-311+g(3df,3pd)': [0.9265, 0.8864, 0.8938, '10.1021/jp073974n'],
        '6-311+g(3df,p)': [0.9262, 0.8878, 0.8953, '10.1021/jp073974n'],
        '6-311+g(d)': [0.9222, 0.8957, 0.903, '10.1021/jp073974n'],
        '6-311+g(d,p)': [0.9255, 0.8967, 0.9041, '10.1021/jp073974n'],
        '6-311g(d)': [0.9214, 0.8947, 0.9018, '10.1021/jp073974n'],
        '6-311g(d,p)': [0.9247, 0.8951, 0.9021, '10.1021/jp073974n'],
        '6-311g(df,p)': [0.9247, 0.8908, 0.8981, '10.1021/jp960976r'],
        '6-31g(2df,p)': [0.9222, 0.8878, 0.8949, '10.1021/jp073974n'],
        '6-31g(d)': [0.9135, 0.8902, 0.8978, '10.1021/jp073974n'],
        '6-31g(d,p)': [0.9181, 0.8913, 0.8991, '10.1021/jp073974n'],
        'tzv2p': [0.9256, 0.8979, 0.9054, '10.1021/jp073974n'],
        'aug-cc-pvdz': [0.9232, 1.1169, 1.1254, '10.1021/jp048233q'],
        'aug-cc-pvqz': [0.9211, 1.1272, 1.1352, '10.1021/jp048233q'],
        'aug-cc-pvtz': [0.9216, 1.1259, 1.1342, '10.1021/jp048233q'],
        'cc-pvdz': [0.9204, 1.1254, 1.137, '10.1021/jp048233q'],
        'cc-pvqz': [0.9202, 1.1284, 1.137, '10.1021/jp048233q'],
        'cc-pvtz': [0.9213, 1.1278, 1.1363, '10.1021/jp048233q']},
 'm05': {'6-31+g(d,p)': [0.9787, 0.9761, 0.9771, '10.1021/jp073974n'],
         '6-311+g(2df,p)': [0.9851, 0.9571, 0.9548, '10.1021/jp073974n'],
         '6-311+g(d,p)': [0.9841, 0.9752, 0.9757, '10.1021/jp073974n'],
         '6-31g(2df,p)': [0.9809, 0.9558, 0.9555, '10.1021/jp073974n'],
         '6-31g(d)': [0.9736, 0.9712, 0.9715, '10.1021/jp073974n'],
         'pc-0': [0.977, 0.9983, 0.9772, '10.1002/jcc.23073'],
         'pc-1': [0.9751, 1.047, 1.056, '10.1002/jcc.23073'],
         'pc-2': [0.9771, 1.046, 1.05, '10.1002/jcc.23073'],
         'pc-3': [0.9806, 1.042, 1.043, '10.1002/jcc.23073'],
         'pc-4': [0.9848, 0.9916, 0.9916, '10.1002/jcc.23073']},
 'm052x': {'6-31+g(d,p)': [0.9631, 0.9504, 0.9445, '10.1021/jp073974n'],
            '6-311+g(2df,p)': [0.9663, 0.9297, 0.9206, '10.1021/jp073974n'],
            '6-311+g(d,p)': [0.9658, 0.9467, 0.9398, '10.1021/jp073974n'],
            '6-31g(2df,p)': [0.9657, 0.9321, 0.9239, '10.1021/jp073974n'],
            '6-31g(d)': [0.958, 0.947, 0.9425, '10.1021/jp073974n'],
            'pc-0': [0.9635, 1.006, 0.97, '10.1002/jcc.23073'],
            'pc-1': [0.9656, 1.047, 1.043, '10.1002/jcc.23073'],
            'pc-2': [0.9661, 1.056, 1.052, '10.1002/jcc.23073'],
            'pc-3': [0.9676, 1.051, 1.041, '10.1002/jcc.23073'],
            'pc-4': [0.967, 1.036, 1.036, '10.1002/jcc.23073']},
 'm06': {'def2-tzvp': [0.9820, 1., 1., 'truhlar'],
        'def2-tzvpp': [0.9790, 1., 1., 'truhlar'],
        'pc-0': [0.986, 0.9945, 0.9744, '10.1002/jcc.23073'],
        'pc-1': [0.9863, 1.042, 1.052, '10.1002/jcc.23073'],
        'pc-2': [0.9849, 1.05, 1.056, '10.1002/jcc.23073'],
        'pc-3': [0.9829, 1.044, 1.045, '10.1002/jcc.23073'],
        'pc-4': [0.9853, 1.007, 1.007, '10.1002/jcc.23073']},
 'm062x': {'def2-tzvp': [0.9710, 1., 1., 'truhlar'],
        'def2-tzvpp': [0.9700, 1., 1., 'truhlar'],
        'def2-tqzvp': [0.9700, 1., 1., 'truhlar'],
        '6-31+g(d,p)': [0.9670, 1., 1., 'truhlar'],
        '6-311+g(d,p)': [0.9700, 1., 1., 'truhlar'],
        'pc-0': [0.9737, 0.9969, 0.9691, '10.1002/jcc.23073'],
        'pc-1': [0.9728, 1.05, 1.052, '10.1002/jcc.23073'],
        'pc-2': [0.9728, 1.059, 1.061, '10.1002/jcc.23073'],
        'pc-3': [0.9734, 1.054, 1.054, '10.1002/jcc.23073'],
        'pc-4': [0.9733, 1.055, 1.054, '10.1002/jcc.23073']},
 'mp2': {'6-31+g(d,p)': [0.9657, 1.0198, 1.0334, '10.1021/jp073974n'],
         '6-311+g(2df,p)': [0.9777, 0.9826, 0.9908, '10.1021/jp073974n'],
         '6-311+g(d,p)': [0.9768, 1.0071, 1.0158, '10.1021/jp073974n'],
         '6-31g(2df,p)': [0.9678, 0.9823, 0.9858, '10.1021/jp073974n'],
         '6-31g(d)': [0.967, 1.0059, 1.0178, '10.1021/jp073974n'],
         'aug-cc-pvdz': [0.9675, 0.9473, 0.9049, '10.1021/jp048233q'],
         'aug-cc-pvqz': [0.9959, 0.9658, 0.9264, '10.1021/jp048233q'],
         'aug-cc-pvtz': [0.983, 0.9901, 0.9467, '10.1021/jp048233q'],
         'cc-pvdz': [0.9784, 0.9796, 0.9379, '10.1021/jp048233q'],
         'cc-pvqz': [0.9943, 0.9957, 0.9575, '10.1021/jp048233q'],
         'cc-pvtz': [0.9832, 0.9907, 0.9477, '10.1021/jp048233q']},
 'o3lyp': {'6-31+g(d,p)': [0.9867, 1.0098, 1.0139, '10.1021/jp073974n'],
           '6-311+g(2df,p)': [0.9918, 0.9965, 1.001, '10.1021/jp073974n'],
           '6-311+g(d,p)': [0.9904, 1.0123, 1.0186, '10.1021/jp073974n'],
           '6-31g(2df,p)': [0.9872, 0.995, 0.9989, '10.1021/jp073974n'],
           '6-31g(d)': [0.9826, 1.0048, 1.0075, '10.1021/jp073974n']},
 'pbe0': {'6-31+g(d,p)': [0.9771, 0.9827, 0.9875, '10.1021/jp073974n'],
          '6-311+g(2df,p)': [0.9824, 0.972, 0.9776, '10.1021/jp073974n'],
          '6-311+g(d,p)': [0.9812, 0.9856, 0.9915, '10.1021/jp073974n'],
          '6-31g(2df,p)': [0.9779, 0.9703, 0.9755, '10.1021/jp073974n'],
          '6-31g(d)': [0.9726, 0.9777, 0.9814, '10.1021/jp073974n']},
 'pbe1pbe': {'pc-0': [0.9793, 0.9911, 0.9693, '10.1002/jcc.23073'],
             'pc-1': [0.9784, 1.045, 1.049, '10.1002/jcc.23073'],
             'pc-2': [0.9812, 1.043, 1.043, '10.1002/jcc.23073'],
             'pc-3': [0.9822, 1.042, 1.042, '10.1002/jcc.23073'],
             'pc-4': [0.9823, 1.042, 1.042, '10.1002/jcc.23073']},
 'pbepbe': {'pc-0': [1.015, 0.9257, 0.9036, '10.1002/jcc.23073'],
            'pc-1': [1.014, 0.9588, 0.9579, '10.1002/jcc.23073'],
            'pc-2': [1.015, 0.9577, 0.9539, '10.1002/jcc.23073'],
            'pc-3': [1.015, 0.9571, 0.9529, '10.1002/jcc.23073'],
            'pc-4': [1.016, 0.9573, 0.9532, '10.1002/jcc.23073']},
 'qcisd': {'6-31+g(d,p)': [0.9703, 1.01, 1.0202, '10.1021/jp073974n'],
           '6-311+g(d,p)': [0.9812, 0.997, 1.0049, '10.1021/jp073974n'],
           '6-31g(d)': [0.9777, 0.9986, 1.0058, '10.1021/jp073974n']},
 'qcisd(t)': {'6-31+g(d,p)': [0.9786, 1.0419, 1.0583, '10.1021/jp073974n'],
              '6-311+g(d,p)': [0.9907, 1.0274, 1.0382, '10.1021/jp073974n'],
              '6-31g(d)': [0.9859, 1.0273, 1.0405, '10.1021/jp073974n']},
 'tpss': {'6-31+g(d,p)': [0.9957, 1.0424, 1.0481, '10.1021/jp073974n'],
          '6-311+g(2df,p)': [1.0007, 1.0263, 1.0308, '10.1021/jp073974n'],
          '6-311+g(d,p)': [0.9999, 1.0463, 1.0554, '10.1021/jp073974n'],
          '6-31g(2df,p)': [0.9971, 1.0264, 1.0308, '10.1021/jp073974n'],
          '6-31g(d)': [0.9925, 1.0382, 1.0425, '10.1021/jp073974n']},
 'vsxc': {'6-31+g(d,p)': [0.9904, 1.0393, 1.0413, '10.1021/jp073974n'],
          '6-311+g(2df,p)': [0.9947, 1.0174, 1.0151, '10.1021/jp073974n'],
          '6-311+g(d,p)': [0.9937, 1.0406, 1.0435, '10.1021/jp073974n'],
          '6-31g(2df,p)': [0.9896, 1.0164, 1.0156, '10.1021/jp073974n'],
          '6-31g(d)': [0.9877, 1.0359, 1.0374, '10.1021/jp073974n']}
}

def MethodPraser(method:str):
    # not implemented
    # converting CCSD-T to CCSD(T) etc.
    return method

class GaussianOutput():
    def __init__(self, filename:str) -> None:
        print(f"\tNow reading {filename}")
        data = cclib.io.ccopen(filename).parse()

        # Check if the Gaussian terminate nomally
        if not data.metadata["success"]:
            print("Error termination! Gaussian did not end normally.")
            return None

        self.basis_set = data.metadata["basis_set"]

        if "methods" in data.metadata.keys():
            self.method = data.metadata["methods"][-1]
            if self.method == "DFT":
                self.method = data.metadata["functional"]
            else:
                self.method = MethodPraser(self.method)

        # get energy
        if "ccenergies" in dir(data):
            self.energy = data.ccenergies[-1] / 27.2114 # convert eV to Hartree
        elif "mpenergies" in dir(data):
            self.energy = data.mpenergies[-1] / 27.2114
        elif "scfenergies" in dir(data):
            self.energy = data.scfenergies[-1] / 27.2114
        else:
            raise AttributeError()

        print(f"\tComputation was conduct in {self.method}/{self.basis_set} level")
        print(f"\tEnergy = {self.energy}")   

class MoleculeSummarize():
    def __init__(self, name) -> None:
        self.name = name

    def __Shermo_parser(self, Shermo_output:str):
        pattern_cv = r"Total CV:\s+([\.0-9]+)\s+J/mol/K\s+([\.0-9]+)\s+cal/mol/K"
        pattern_cp = r"Total CP:\s+([\.0-9]+)\s+J/mol/K\s+([\.0-9]+)\s+cal/mol/K"
        pattern_s = r"Total S:\s+([\.0-9]+)\s+J/mol/K\s+([\.0-9]+)\s+cal/mol/K"
        pattern_zpe = r"Zero point energy \(ZPE\):\s+([\.0-9]+)\s+kJ/mol\s+([\.0-9]+)\s+kcal/mol\s+([\.0-9]+)\s+a\.u\."
        pattern_c_u = r"Thermal correction to U:\s+([\.0-9]+)\s+kJ/mol\s+([\.0-9]+)\s+kcal/mol\s+([\.0-9]+)\s+a\.u\."
        pattern_c_h = r"Thermal correction to H:\s+([\.0-9]+)\s+kJ/mol\s+([\.0-9]+)\s+kcal/mol\s+([\.0-9]+)\s+a\.u\."
        pattern_c_g = r"Thermal correction to G:\s+([\.0-9]+)\s+kJ/mol\s+([\.0-9]+)\s+kcal/mol\s+([\.0-9]+)\s+a\.u\."

        self.cv = float(re.search(pattern_cv, Shermo_output).group(1)) # group 1 in J/mol/K; group 2 in cal/mol/K
        self.cp = float(re.search(pattern_cp, Shermo_output).group(1)) # group 1 in J/mol/K; group 2 in cal/mol/K
        self.s = float(re.search(pattern_s, Shermo_output).group(1)) # group 1 in J/mol/K; group 2 in cal/mol/K
        self.zpe = float(re.search(pattern_zpe, Shermo_output).group(3)) # group 1 in kJ/mol; group 2 in kcal/mol; group 3 in a.u.
        self.c_u = float(re.search(pattern_c_u, Shermo_output).group(3)) # group 1 in kJ/mol; group 2 in kcal/mol; group 3 in a.u.
        self.c_h = float(re.search(pattern_c_h, Shermo_output).group(3)) # group 1 in kJ/mol; group 2 in kcal/mol; group 3 in a.u.
        self.c_g = float(re.search(pattern_c_g, Shermo_output).group(3)) # group 1 in kJ/mol; group 2 in kcal/mol; group 3 in a.u.
        
        self.u = self.energy + self.c_u
        self.h = self.energy + self.c_h
        self.g = self.energy + self.c_g

    def getEnergy(self, log:str):
        data = GaussianOutput(log)
        self.energy_log = log
        self.energy = data.energy
        self.method_energy = data.method + "/" + data.basis_set

    def getFreq(self, log:str):
        data = GaussianOutput(log)
        self.freq_log = log
        method = data.method.lower()
        basis_set = data.basis_set.lower()
        self.method_freq = data.method + "/" + data.basis_set

        try:
            # actually alternative scaler can be use (no diffuse or no polarization )
            self.scaler = ThermoScaler[method][basis_set]
            print("Using scaler", self.scaler)
        except KeyError:
            print("Did not find suitable scaler")
            self.scaler = [1., 1., 1., "No scaler"]

    def runShermo(self, output_dir) -> int:
        Shermo_params = [shermo_path, self.freq_log,
                      "-E", str(self.energy),
                      "-sclZPE", str(self.scaler[0]),
                      "-sclheat", str(self.scaler[1]),
                      "-sclS", str(self.scaler[2])
                      ]
        result = subprocess.run(Shermo_params, capture_output=True, text=True)

        output_file = os.path.join(output_dir, f"{self.name}.txt")
        if result.returncode == 0:
            print(f"Shermo completed successfully for molecule {self.name}.")
            self.__Shermo_parser(result.stdout)
            with open(output_file, 'w+') as f:
                f.write(result.stdout)
            return 0
        else:
            print(f"Shermo failed for molecule {self.name}.")
            with open(output_file, 'w+') as f:
                f.write(result.stderr)
            return 1



def ProjectSummarize(cwd=os.getcwd()):
    def getLog(files):
    # get .log files during os.walk
        log_files_list = []
        for file in files:
            if file.endswith('.log'):
                full_file_dir = os.path.join(root, file)
                print("Found:", full_file_dir)
                log_files_list.append(full_file_dir)
        return log_files_list
    
    # Looking for energy, opt, and freq directory
    for root, dirs, files in os.walk(cwd):
        if os.path.basename(root) == DIR_OPT:
            print("\nNow looking for optimizations")
            opt_files = getLog(files)
        elif os.path.basename(root) == DIR_FREQ:
            print("\nNow looking for frequency")
            freq_files = getLog(files)
        elif os.path.basename(root) == DIR_ENERGY:
            print("\nNow looking for energy")
            energy_files = getLog(files)

    # Make a dir for Shermo output
    shermo_output_dir = os.path.join(cwd, "ShermoOutput")
    if not os.path.exists(shermo_output_dir):
        os.makedirs(shermo_output_dir)

    # Run Shermo for each molecule found in energy dir
    mols = []
    for energy_file in energy_files:
        name = os.path.basename(energy_file).split("_")[0]
        print("\nNow reading molecule", name)
        
        mol = MoleculeSummarize(name)
        print("Getting precise single point energy calculation:")
        mol.getEnergy(energy_file)

        print("Getting frequency calculation")
        # try looking for freq computation in freq dir
        found_in_freq_dir = False
        for freq_file in freq_files:
            if os.path.basename(freq_file).split("_")[0] == name:
                mol.getFreq(freq_file)
                found_in_freq_dir = True
                if mol.runShermo(output_dir=shermo_output_dir) == 0:
                    mols.append(mol)

        # some tasks run freq calculation during opt
        if not found_in_freq_dir:
            for opt_file in opt_files:
                if os.path.basename(opt_file).split(".")[0] == name:
                    mol.getFreq(opt_file)
                    if mol.runShermo(output_dir=shermo_output_dir) == 0:
                        mols.append(mol)

    with open(os.path.join(shermo_output_dir, "summay.csv"), "w+", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Name",
                         "Method for freq",
                         "Method for energy",
                         "Electronic energy (a.u.)",
                         "Cv (J/mol/K)",
                         "Cp (J/mol/K)",
                         "Entropy (J/mol/K)",
                         "ZPE (a.u.)",
                         "Thermal correction to U (a.u.)",
                         "Thermal correction to H (a.u.)",
                         "Thermal correction to G (a.u.)",
                         "EE + ZPE (a.u.)",
                         "EE + correction to U (a.u.)",
                         "EE + correction to H (a.u.)",
                         "EE + correction to G (a.u.)",
                         "Thermo scaler"])
        for mol in mols:
            writer.writerow([mol.name, mol.method_freq, mol.method_energy,
                             mol.energy, mol.cv, mol.cp, mol.s,
                             mol.zpe, mol.c_u, mol.c_h, mol.c_g,
                             mol.zpe+mol.energy, mol.u, mol.h, mol.g,
                             mol.scaler])

if __name__=="__main__":
    ProjectSummarize()