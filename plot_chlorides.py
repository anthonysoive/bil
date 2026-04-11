import numpy as np
import matplotlib.pyplot as plt
import re
import os

fileName = "ChloricemDLM.t12"

with open(fileName, 'r') as f:
    lines = f.readlines()
    
header_line = ""
for line in lines:
    if line.startswith("# Coordinates(1)"):
        header_line = line[2:].strip()
        break

pattern = re.compile(r"([^()]+)\(\d+\)")
matches = [m.strip() for m in pattern.findall(header_line)]

idx_x = matches.index('Coordinates')
idx_tot = matches.index('n_Cl')
idx_ads = matches.index('adsorbed chloride')
idx_fs = matches.index("n_Friedel's salt")

data = np.loadtxt(fileName, comments='#')
# depth in cm
depth = data[:, idx_x] * 100 # from decimeter to mm? No, 1 dm = 10 cm. If we want cm, * 10. If mm, * 100.
depth = data[:, idx_x] * 10 # cm

tot_cl = data[:, idx_tot]
ads_cl = data[:, idx_ads]
fs_cl = 2 * data[:, idx_fs]
free_cl = tot_cl - ads_cl - fs_cl

plt.figure(figsize=(10,6))
# Tracer les aires cumulées avec fill_between
plt.fill_between(depth, 0, free_cl, color='royalblue', alpha=0.6, label='Chlorures libres')
plt.fill_between(depth, free_cl, free_cl + ads_cl, color='mediumseagreen', alpha=0.6, label='Chlorures adsorbés (C-S-H)')
plt.fill_between(depth, free_cl + ads_cl, free_cl + ads_cl + fs_cl, color='indianred', alpha=0.6, label='Chlorures fixés (Sel de Friedel)')

# Ligne du total pour bien marquer la limite
plt.plot(depth, tot_cl, 'k-', linewidth=1.5, label='Chlorures totaux (courbe)')

plt.xlabel('Profondeur (cm)')
plt.ylabel('Concentration en chlorures (mol/L de matériau)')
plt.title('Profil de pénétration des chlorures dans le béton (1 an)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('ChloricemDLM_results.png', dpi=300, bbox_inches='tight')
print("Plot saved to ChloricemDLM_results.png")
