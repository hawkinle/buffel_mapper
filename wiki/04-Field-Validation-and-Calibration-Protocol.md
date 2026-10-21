# 04. Field Validation & Calibration Protocol

**Navigation:** [Home](Home.md) | [01. Methodology](01-Ecology-and-Remote-Sensing-Methodology.md) | [02. Architecture & CLI](02-Pipeline-Architecture-and-CLI-Reference.md) | [03. Compilation](03-Compilation-and-Deployment-Guide.md) | **04. Field Validation**

---

## 1. Field Schema & Classification Rules

To validate and calibrate the 3-Tier remote sensing model in the paddock, record observations at each GPS validation point within a **5-meter radius** centered on the coordinates:

### A. Weed Cover Class (`Observed_Weed_Cover_Class_0_4`)
| Class | Description | Foliage Cover % | Operational Tier Equivalent |
| :---: | :--- | :---: | :--- |
| **0** | **Absent** | $0\%$ | Background / Clean (`0`) |
| **1** | **Rare / Isolated** | $< 5\%$ | Tier 1 — Eradication (`1`) |
| **2** | **Sparse / Scattered** | $5\% - 20\%$ | Tier 1 — Eradication (`1`) |
| **3** | **Moderate** | $20\% - 50\%$ | Tier 2 — Containment (`2`) |
| **4** | **Dense Core** | $> 50\%$ | Tier 2 — Containment (`2`) |

### B. Overstory Canopy Class (`Observed_Overstory_Canopy_Class_0_3`)
| Class | Description | Canopy Cover % | Refugia Delineation |
| :---: | :--- | :---: | :--- |
| **0** | **Absent** | $0\%$ | Open grassland / plain |
| **1** | **Very Sparse** | $< 10\%$ | Scattered isolated trees |
| **2** | **Open Canopy** | $10\% - 30\%$ | Open woodland (Eucalyptus, Mulga, Rosewood) |
| **3** | **Dense Canopy** | $> 30\%$ | Validated Tier 3 Native Refugia (`3`) |

---

## 2. Field Survey Procedure

1. **Load Waypoints:** Load the validation points (`.gpx`, `.kml`, or `.csv`) into a mobile GIS app (Avenza Maps, QField, or handheld GPS).
2. **Navigate & Photograph:** Stand at the coordinate center and capture 4 cardinal photos (North, East, South, West).
3. **Record Classes:** Populate `Observed_Weed_Cover_Class_0_4` and `Observed_Overstory_Canopy_Class_0_3` in the validation CSV.
4. **Run Calibration & Accuracy Audit:**
   ```bash
   ./dist/buffel-mapper calibrate \
     --csv-path /path/to/field_validation_points.csv \
     --project-dir ./projects/Keewong_Run \
     --optimize
   ```
   The calibration script outputs an overall accuracy score, confusion matrix, and optimized spectral thresholds tailored to your ground observations.
