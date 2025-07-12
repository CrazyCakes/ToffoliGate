This is the code and the raw data of the paper: Generalized Toffoli Gates with Customizable Single-Step Multiple-Qubit Control by 
Chung-Kai Wu, Dah-Wei Chiou, and Jie-Hong Roland Jiang.
Graduate Institute of Electronics Engineering, National Taiwan University, Taipei, Taiwan

This repository contains:

| File / Folder   | Description                                                                                                   |
|-----------------|---------------------------------------------------------------------------------------------------------------|
| `Simulation.py` | QuTiP‑based script that evaluates average gate fidelities for multi‑qubit Toffoli‑type operations under various hardware and noise settings. |
| `RawData.xlsx`  | Spreadsheet with all simulated and experimental data used for the figures and tables in the paper.           |
| `README.md`     | Project documentation (this file).                                                                           |

---


### Clone and set up the environment

```bash
git clone https://github.com/CrazyCakes/ToffoliGate.git
cd <repo‑name>

# (optional) create a virtual environment
python -m venv .venv
source .venv/bin/activate      # Linux / macOS
# .venv\Scripts\activate.bat  # Windows (PowerShell)

pip install numpy qutip
```


###  Run the simulation with default settings

```bash
python3 Simulation.py
```
