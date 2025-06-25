# Import necessary modules
from rs2.modeler.RS2Modeler import RS2Modeler
from rs2.modeler.properties.PropertyEnums import *
from rs2.interpreter.RS2Interpreter import RS2Interpreter
from rs2.interpreter.InterpreterEnums import *
import os
import math
import time
from datetime import datetime

# --- Inputs ---
MODEL_FILE = r"Absolute path to the RS2 input model file"
BASE_OUTPUT_DIR = r"Absolute path to the output directory where results will be saved"

# Material parameters
GSI = 65
E_i = 18800  # MPa
m_i = 10.75
D = 0.8
UCS = 67.48  # MPa
Poisson = 0.15
Porosity = 0.1
UnitWeight = 0.0285  # MN/m^3
Name_of_material = "Dunderland_Schist"

# Target displacement and tolerance
target_change_mm = -15.0
tolerance_mm = 0.1 * abs(target_change_mm)
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

# Number of iterations and ports used
max_iterations = 50
modeler_port = 60502
interpreter_port = 60503

# --- End of inputs ---

start_time = time.time()

# Calculate original distance between points
original_distance = abs(point2[0] - point1[0]) # 4.3 meters

# --- Helper Functions ---
def calculate_hoek_brown_parameters(GSI, m_i, D):
    mb = m_i * math.exp((GSI - 100) / (28 - 14 * D))
    s = math.exp((GSI - 100) / (9 - 3 * D))
    a = 0.5 + (1/6) * (math.exp(-GSI / 15) - math.exp(-20 / 3))
    return mb, s, a

def calculate_Em(GSI, E_i, D):
    numerator = (1 - D/2)
    denominator = 1 + math.exp((60 + 15*D - GSI) / 11)
    E_m = E_i * (0.02 + numerator / denominator)
    return E_m

def run_model_and_get_change(GSI, E_i, D, modeler):
    mb, s, a = calculate_hoek_brown_parameters(GSI, m_i, D)
    E_m = calculate_Em(GSI, E_i, D)

    print("Opening model...", end=" ")
    model = modeler.openFile(MODEL_FILE)

    # Set material properties
    print("Setting material properties...", end=" ")
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

    model_path = os.path.join(OUTPUT_DIR, "U_stoll_Elastic.fez")
    if os.path.exists(model_path):
        os.remove(model_path)

    print("Saving model...", end=" ")
    model.saveAs(model_path)
    print("Computing model...", end=" ")
    model.compute()
    print("Analysis done. Getting displacements...", end=" ")

    # Get displacement results using Pythagoras' theorem
    RS2Interpreter.startApplication(port=interpreter_port)
    interpreter = RS2Interpreter(port=interpreter_port)
    results = interpreter.openFile(model_path)

    results.AddMaterialQuery(points=[point1, point2])
    results.SetResultType(ExportResultType.SOLID_DISPLACEMENT_HORIZONTAL_DISPLACEMENT)
    ux1, ux2 = [p.GetValue() for p in results.GetMaterialQueryResults()[0].GetAllValues()]
    results.SetResultType(ExportResultType.SOLID_DISPLACEMENT_VERTICAL_DISPLACEMENT)
    uy1, uy2 = [p.GetValue() for p in results.GetMaterialQueryResults()[0].GetAllValues()]

    delta_x = (point2[0] - point1[0]) + (ux2 - ux1)
    delta_y = (point2[1] - point1[1]) + (uy2 - uy1)
    new_distance = math.sqrt(delta_x**2 + delta_y**2)
    change_mm = (new_distance - original_distance) * 1000

    model.close()
    try:
        interpreter.closeProgram()
        print("Interpreter closed.")
    except Exception as e:
        print(f"Failed to close Interpreter: {e}")

    return change_mm, mb, s, a, E_m

# Start RS2 Modeler
try:
    print("Starting RS2 Modeler...")
    RS2Modeler.startApplication(port=modeler_port)
    modeler = RS2Modeler(port=modeler_port)
except Exception as e:
    print(f"Failed to start RS2 Modeler: {e}")
    raise

# Storage lists
valid_results = []
change_history = []
gsi_history = []
ei_history = []
within_range = False
entry_direction = None

# Simulation loop
try:
    gsi_iter = GSI
    ei_iter = E_i

    for iteration in range(max_iterations):
        print(f"\nIteration {iteration + 1}")
        change_mm, mb, s, a, E_m = run_model_and_get_change(gsi_iter, ei_iter, D, modeler)
        change_history.append(change_mm)
        gsi_history.append(gsi_iter)
        ei_history.append(ei_iter)

        print(f"Used in iteration {iteration + 1} | GSI: {gsi_iter}, E_i: {ei_iter}, Change: {change_mm:.2f} mm")

        if min_change_mm <= change_mm <= max_change_mm:
            if not within_range and iteration > 0:
                prev = change_history[-2]
                entry_direction = "too_negative" if prev < min_change_mm else "too_positive"
                print(f"Entered range from {entry_direction.replace('_', ' ')} side")
            within_range = True
            valid_results.append((gsi_iter, ei_iter, change_mm, mb, s, a, E_m))
        else:
            if within_range:
                print("Outside target range. Stopping.")
                break

        dist = abs(change_mm - (min_change_mm if change_mm < min_change_mm else max_change_mm))
        if not within_range:
            step_GSI, step_E_i = (6, 600) if dist > 3 else (2, 200) if dist > 1 else (1, 100)
        else:
            step_GSI, step_E_i = (3, 300) if dist > 1 else (1, 100)

        if (iteration + 1) % 2 == 1:
            gsi_iter += step_GSI if (not within_range and change_mm < min_change_mm) or (within_range and entry_direction == "too_negative") else -step_GSI
            gsi_iter = max(min(gsi_iter, GSI_max), GSI_min)
        else:
            ei_iter += step_E_i if (not within_range and change_mm < min_change_mm) or (within_range and entry_direction == "too_negative") else -step_E_i
            ei_iter = max(min(ei_iter, E_i_max), E_i_min)

        print(f" -> Next step: GSI={gsi_iter}, E_i={ei_iter}")
        print("-----")

    print("\nSimulation finished.")
    print(f"Total iterations: {len(change_history)}")

finally:
    try:
        modeler.closeProgram()
        print("RS2 Modeler closed successfully.")
    except Exception as e:
        print(f"Failed to close RS2 Modeler: {e}")

    total_time = time.time() - start_time
    minutes = int(total_time // 60)
    seconds = int(total_time % 60)

    output_lines = []
    output_lines.append("Input Parameters")
    output_lines.append("----------------")
    output_lines.append(f"Material Parameters: GSI={GSI}, E_i={E_i} MPa, m_i={m_i}, D={D}, UCS={UCS} MPa, Poisson={Poisson}, Porosity={Porosity}, UnitWeight={UnitWeight} MN/m^3, Material={Name_of_material}")
    output_lines.append(f"Target Displacement: {target_change_mm} mm, Tolerance: 10.0%, Min Change: {min_change_mm} mm, Max Change: {max_change_mm} mm")
    output_lines.append(f"Points for Measuring Convergence: point1={point1}, point2={point2}")
    output_lines.append(f"Parameter Ranges: GSI_min={GSI_min}, GSI_max={GSI_max}, E_i_min={E_i_min:.1f}, E_i_max={E_i_max:.1f}")
    output_lines.append(f"Simulation Settings: Iterations={len(change_history)}, modeler_port={modeler_port}, interpreter_port={interpreter_port}")
    output_lines.append("")
    output_lines.append("Optimization finished.")
    output_lines.append(f"Time taken: {minutes} min and {seconds} sec\n")

    output_lines.append("All iterations:")
    output_lines.append("--------------------------------------------------")
    for i in range(len(change_history)):
        output_lines.append(f"Iteration {i + 1}: GSI={gsi_history[i]:.3f}, E_i={ei_history[i]:.3f}, Displacement={change_history[i]:.3f} mm")

    if valid_results:
        output_lines.append("\nBest combinations within tolerance:")
        output_lines.append("--------------------------------------------------")
        best = sorted(valid_results, key=lambda x: abs(x[2] - target_change_mm))[:5]
        for gsi, ei, disp, mb, s, a, em in best:
            output_lines.append(f"GSI={gsi:.3f}, E_i={ei:.3f}, Displacement={disp:.3f} mm, mb={mb:.3f}, s={s:.3f}, a={a:.3f}, E_m={em:.3f}")

        output_lines.append("\nSpecific combinations:")
        output_lines.append("--------------------------------------------------")
        lowest_gsi_idx = valid_results.index(min(valid_results, key=lambda x: x[0]))
        highest_gsi_idx = valid_results.index(max(valid_results, key=lambda x: x[0]))
        lowest_ei_idx = valid_results.index(min(valid_results, key=lambda x: x[1]))
        highest_ei_idx = valid_results.index(max(valid_results, key=lambda x: x[1]))

        for label, idx in [
            ("Lowest GSI", lowest_gsi_idx),
            ("Highest GSI", highest_gsi_idx),
            ("Lowest E_i", lowest_ei_idx),
            ("Highest E_i", highest_ei_idx),
        ]:
            gsi, ei, disp, *_ = valid_results[idx]
            output_lines.append(f"{label}: GSI={gsi:.3f}, E_i={ei:.3f}, Displacement={disp:.3f} mm")

    filename = f"Elastic_Hardcoded_{datetime.now().strftime('%Y%m%d')}_Iterations{len(change_history)}.txt"
    result_path = os.path.join(OUTPUT_DIR, filename)
    with open(result_path, "w") as f:
        for line in output_lines:
            f.write(line + "\n")
    print(f"Summary saved to: {filename}")

