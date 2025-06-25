# Import necessary modules
from rs2.modeler.RS2Modeler import RS2Modeler 
from rs2.modeler.properties.PropertyEnums import * 
from rs2.interpreter.RS2Interpreter import RS2Interpreter 
from rs2.interpreter.InterpreterEnums import * 

# General-purpose standard library imports
import os 
import math 
import time 
import optuna 
from optuna.samplers import TPESampler 
import optuna.visualization as vis 

# Date and time handling
from datetime import datetime 

# --- Inputs ---
MODEL_FILE = r"Absolute path to the RS2 input model file (.fez)"
BASE_OUTPUT_DIR = r"Absolute path to the output directory where results will be saved"

# Material parameters
GSI = 65
E_i = 18800  # (MPa)
m_i = 10.75
D = 0.8
UCS = 67.48  # (MPa)
Poisson = 0.15
Porosity = 0.1
UnitWeight = 0.0285  # (MN/m^3)
Name_of_material = "Dunderland_Schist"

# Target displacement and tolerance
target_change_mm = -15.0  # Target change in distance (mm)
tolerance_percent = 10.0  # Tolerance as a percentage of target_change_mm
tolerance_mm = (tolerance_percent / 100) * abs(target_change_mm)
min_change_mm = target_change_mm - tolerance_mm  
max_change_mm = target_change_mm + tolerance_mm 

# Points for measuring convergence
point1 = (22.850, 23.966)
point2 = (27.150, 23.966)

# Parameter ranges
GSI_min = GSI - 20
GSI_max = GSI + 20
E_i_min = E_i * 0.9  
E_i_max = E_i * 1.1  

# Simulation settings
n_trials = 100
modeler_port = 60500
interpreter_port = 60501

# TPESampler parameters (adjustable)
n_startup_trials = 10 # Default is 10
n_ei_candidates = 24 # Default is 24

# --- End of inputs ---

# Create a unique directory for this run
current_date = datetime.now().strftime("%Y%m%d")
RUN_DIR = os.path.join(BASE_OUTPUT_DIR, f"Elastic_Optuna_{current_date}_Trials{n_trials}")
os.makedirs(RUN_DIR, exist_ok=True)

# Define results file within the run directory
RESULTS_FILE = os.path.join(RUN_DIR, f"Optuna_results_{current_date}_{n_trials}.txt")

# Calculate original distance between points
original_distance = abs(point2[0] - point1[0])  # 4.3 meters

# --- Helper Functions ---
def calculate_hoek_brown_parameters(GSI, m_i, D):
    """Calculate Hoek-Brown parameters mb, s, and a based on GSI, m_i, and D."""
    mb = m_i * math.exp((GSI - 100) / (28 - 14 * D))
    s = math.exp((GSI - 100) / (9 - 3 * D))
    a = 0.5 + (1/6) * (math.exp(-GSI / 15) - math.exp(-20 / 3))
    return mb, s, a

def calculate_Em(GSI, E_i, D):
    """Calculate rock mass modulus E_m based on GSI, E_i, and D."""
    numerator = (1 - D/2)
    denominator = 1 + math.exp((60 + 15*D - GSI) / 11)
    E_m = E_i * (0.02 + numerator / denominator)
    return E_m

def run_model_and_get_change(GSI, E_i, D, modeler):
    """Run RS2 model with given parameters and return displacement change in mm."""
    mb, s, a = calculate_hoek_brown_parameters(GSI, m_i, D)
    E_m = calculate_Em(GSI, E_i, D)

    print(f"Trial with GSI={GSI:.3f}, E_i={E_i:.3f} | Opening model...", end=" ")
    model = modeler.openFile(MODEL_FILE)

    # Set material properties
    Dunderland_Schist = model.getAllMaterialProperties()[0]
    Dunderland_Schist.setMaterialName(Name_of_material)
    Dunderland_Schist.InitialConditions.setUnitWeight(UnitWeight)
    Dunderland_Schist.InitialConditions.setPorosityValue(Porosity)
    Dunderland_Schist.Stiffness.Isotropic.setYoungsModulus(E_m)
    Dunderland_Schist.Stiffness.Isotropic.setPoissonsRatio(Poisson)
    Dunderland_Schist.Strength.setFailureCriterion(StrengthCriteriaTypes.GENERALIZED_HOEK_BROWN)
    Dunderland_Schist.Strength.GeneralizedHoekBrown.setMaterialType(MaterialType.ELASTIC)
    Dunderland_Schist.Strength.GeneralizedHoekBrown.setCompressiveStrength(UCS)
    Dunderland_Schist.Strength.GeneralizedHoekBrown.setMbParameter(mb)
    Dunderland_Schist.Strength.GeneralizedHoekBrown.setSParameter(s)
    Dunderland_Schist.Strength.GeneralizedHoekBrown.setAParameter(a)

    model.compute()
    print("Computed.", end=" ")

    # Get displacement results using Pythagoras' theorem
    RS2Interpreter.startApplication(port=interpreter_port)
    interpreter = RS2Interpreter(port=interpreter_port)
    model_results = interpreter.openFile(MODEL_FILE)

    points_making_line = [point1, point2]
    lineID = model_results.AddMaterialQuery(points=points_making_line)

    model_results.SetResultType(ExportResultType.SOLID_DISPLACEMENT_HORIZONTAL_DISPLACEMENT)
    query_results_x = model_results.GetMaterialQueryResults()[0]
    query_points_x = query_results_x.GetAllValues()
    ux1 = query_points_x[0].GetValue()
    ux2 = query_points_x[1].GetValue()

    model_results.SetResultType(ExportResultType.SOLID_DISPLACEMENT_VERTICAL_DISPLACEMENT)
    query_results_y = model_results.GetMaterialQueryResults()[0]
    query_points_y = query_results_y.GetAllValues()
    uy1 = query_points_y[0].GetValue()
    uy2 = query_points_y[1].GetValue()

    delta_x = (point2[0] - point1[0]) + (ux2 - ux1)
    delta_y = (point2[1] - point1[1]) + (uy2 - uy1)
    new_distance = math.sqrt(delta_x**2 + delta_y**2)
    change_in_distance_mm = (new_distance - original_distance) * 1000

    print(f" Displacement: {change_in_distance_mm:.3f} mm")

    model.close()
    interpreter.closeProgram()

    return change_in_distance_mm, mb, s, a, E_m

# --- Main Execution ---
print("Starting RS2 Modeler...")
RS2Modeler.startApplication(port=modeler_port)
modeler = RS2Modeler(port=modeler_port)

def objective(trial):
    """Objective function for Optuna to minimize the difference from target displacement."""
    GSI = trial.suggest_float("GSI", GSI_min, GSI_max)
    E_i = trial.suggest_float("E_i", E_i_min, E_i_max)
    change_mm, mb, s, a, E_m = run_model_and_get_change(GSI, E_i, D, modeler)
    
    trial.set_user_attr("change_mm", change_mm)
    trial.set_user_attr("mb", mb)
    trial.set_user_attr("s", s)
    trial.set_user_attr("a", a)
    trial.set_user_attr("E_m", E_m)
    
    return abs(change_mm - target_change_mm)

# Start timer
start_time = time.time()

# Create Optuna study with multivariate TPESampler and adjusted parameters
sampler = TPESampler(multivariate=True, n_startup_trials=n_startup_trials, n_ei_candidates=n_ei_candidates)
study = optuna.create_study(sampler=sampler, direction="minimize")

study.optimize(objective, n_trials=n_trials)

# --- Evaluate Results and Generate Summary Report  ---
best_trials = []
for trial in study.trials:
    change_mm = trial.user_attrs.get("change_mm", None)
    if change_mm is not None and min_change_mm <= change_mm <= max_change_mm:
        best_trials.append({
            "trial": trial.number,
            "GSI": trial.params["GSI"],
            "E_i": trial.params["E_i"],
            "change_mm": change_mm,
            "mb": trial.user_attrs["mb"],
            "s": trial.user_attrs["s"],
            "a": trial.user_attrs["a"],
            "E_m": trial.user_attrs["E_m"],
            "diff": abs(change_mm - target_change_mm)
        })

best_trials.sort(key=lambda x: x["diff"])

# Find specific combinations
lowest_GSI = min(best_trials, key=lambda x: x["GSI"], default=None)
highest_GSI = max(best_trials, key=lambda x: x["GSI"], default=None)
lowest_E_i = min(best_trials, key=lambda x: x["E_i"], default=None)
highest_E_i = max(best_trials, key=lambda x: x["E_i"], default=None)

# Prepare output
output_lines = []

# Add input parameters to the output
output_lines.append("Input Parameters")
output_lines.append("----------------")
output_lines.append(f"Material Parameters: GSI={GSI}, E_i={E_i} MPa, m_i={m_i}, D={D}, UCS={UCS} MPa, Poisson={Poisson}, Porosity={Porosity}, UnitWeight={UnitWeight} MN/m^3, Material={Name_of_material}")
output_lines.append(f"Target Displacement: {target_change_mm} mm, Tolerance: {tolerance_percent}%, Min Change: {min_change_mm} mm, Max Change: {max_change_mm} mm")
output_lines.append(f"Points for Measuring Convergence: point1={point1}, point2={point2}")
output_lines.append(f"Parameter Ranges: GSI_min={GSI_min}, GSI_max={GSI_max}, E_i_min={E_i_min}, E_i_max={E_i_max}")
output_lines.append(f"Simulation Settings: n_trials={n_trials}, modeler_port={modeler_port}, interpreter_port={interpreter_port}")
output_lines.append(f"TPESampler Parameters: n_startup_trials={n_startup_trials}, n_ei_candidates={n_ei_candidates}")
output_lines.append("")

output_lines.append("Optimization finished.")
total_time = time.time() - start_time
minutes = int(total_time // 60)
seconds = int(total_time % 60)
output_lines.append(f"Time taken: {minutes} min and {seconds} sec")

output_lines.append(f"\nAll trials:")
output_lines.append("--------------------------------------------------")
for trial in study.trials:
    if "change_mm" in trial.user_attrs:
        line = (f"Trial {trial.number}: GSI={trial.params['GSI']:.3f}, E_i={trial.params['E_i']:.3f}, "
                f"Displacement={trial.user_attrs['change_mm']:.3f} mm")
        output_lines.append(line)

output_lines.append(f"\nBest combinations within tolerance ({min_change_mm:.1f} mm to {max_change_mm:.1f} mm):")
output_lines.append("--------------------------------------------------")
for trial in best_trials[:5]: 
    line = (f"Trial {trial['trial']}: GSI={trial['GSI']:.3f}, E_i={trial['E_i']:.3f}, "
            f"Displacement={trial['change_mm']:.3f} mm, Diff={trial['diff']:.3f} mm, "
            f"mb={trial['mb']:.3f}, s={trial['s']:.3f}, a={trial['a']:.3f}, E_m={trial['E_m']:.3f}")
    output_lines.append(line)

if lowest_GSI and highest_GSI and lowest_E_i and highest_E_i:
    output_lines.append("\nSpecific combinations:")
    output_lines.append("--------------------------------------------------")
    output_lines.append(f"Lowest GSI: Trial {lowest_GSI['trial']}, GSI={lowest_GSI['GSI']:.3f}, E_i={lowest_GSI['E_i']:.3f}, "
                        f"Displacement={lowest_GSI['change_mm']:.3f} mm")
    output_lines.append(f"Highest GSI: Trial {highest_GSI['trial']}, GSI={highest_GSI['GSI']:.3f}, E_i={highest_GSI['E_i']:.3f}, "
                        f"Displacement={highest_GSI['change_mm']:.3f} mm")
    output_lines.append(f"Lowest E_i: Trial {lowest_E_i['trial']}, GSI={lowest_E_i['GSI']:.3f}, E_i={lowest_E_i['E_i']:.3f}, "
                        f"Displacement={lowest_E_i['change_mm']:.3f} mm")
    output_lines.append(f"Highest E_i: Trial {highest_E_i['trial']}, GSI={highest_E_i['GSI']:.3f}, E_i={highest_E_i['E_i']:.3f}, "
                        f"Displacement={highest_E_i['change_mm']:.3f} mm")

# Print to console
for line in output_lines:
    print(line)

# Save to file
with open(RESULTS_FILE, "w") as f:
    for line in output_lines:
        f.write(line + "\n")

# --- Visualization ---

print("\nSaving visualizations...")

# 1) Optimization History
fig_hist = vis.plot_optimization_history(study)
fig_hist.write_html(os.path.join(RUN_DIR, "optimization_history.html"))

# 2) Parameter Importances
fig_importance = vis.plot_param_importances(study)
fig_importance.write_html(os.path.join(RUN_DIR, "param_importances.html"))

# 3) Contour Plot (standard)
#fig_contour = vis.plot_contour(study, params=["GSI", "E_i"])
#fig_contour.write_html(os.path.join(RUN_DIR, "contour.html"))

# 3) Contour Plot with finer details and corrected axes
fig_contour = vis.plot_contour(study, params=["E_i", "GSI"])  # E_i on x, GSI on y
fig_contour.update_layout(
    title="Contour Plot - E_i vs GSI",
    xaxis_title="E_i (MPa)",
    yaxis_title="GSI",
    font=dict(size=14),
    plot_bgcolor="white",
    coloraxis=dict(
        colorscale="Rainbow",  
        colorbar=dict(
            title=dict(
                text="Displacement Difference (mm)",
                side="right",
                font=dict(size=14),
            ),
            tickfont=dict(size=12),
            tickmode="array",
            tickvals=[0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20],  
            ticktext=["0", "1", "2", "3", "4", "5", "6", "8", "10", "15", "20"],  
        ),
        cmin=0,  
        cmax=20,  
    ),
)
fig_contour.update_traces(
    selector=dict(type="contour"),
    autocontour=False,
    colorscale="Rainbow",  # Match with coloraxis
    contours=dict(
        coloring="heatmap",
        showlabels=True,
        labelfont=dict(size=12, color="white"),
        start=0,
        end=6,  
        size=0.5,  
    ),
    zmin=0,  
    zmax=20,  
)
fig_contour.write_html(os.path.join(RUN_DIR, "contour.html"))

# 4) Slice Plot
fig_slice = vis.plot_slice(study)
fig_slice.write_html(os.path.join(RUN_DIR, "slice.html"))

# 5) Parallel Coordinate Plot
fig_parallel = vis.plot_parallel_coordinate(study)
fig_parallel.write_html(os.path.join(RUN_DIR, "parallel_coordinate.html"))

# 6) Empirical Distribution Function
fig_edf = vis.plot_edf(study)
fig_edf.write_html(os.path.join(RUN_DIR, "edf.html"))

print("Visualizations saved as HTML in:", RUN_DIR)


# Close RS2 Modeler
modeler.closeProgram()
print("RS2 Modeler closed successfully.")
