import numpy as np
import subprocess
import os
import time

files_to_check = [
    'final_corr1.txt', 'final_corr2.txt', 'final_corr3.txt',
    'final_corr4.txt', 'final_corr5.txt', 'final_corr6.txt'
]

start_time = time.time()
for i in np.linspace(0.03,0.12, 10):
    INPR = i
    np.savetxt('INPR.txt', [INPR], fmt='%.6f')
    result = subprocess.run([r'scripts\betaLoop.bat'], capture_output=True, text=True)
    if result.returncode != 0:
        print("error CorrLoop:", result.stderr)
    else:
        print("Succesfull CorrLoop.")

print("Time taken for the loop: ", time.time() - start_time)