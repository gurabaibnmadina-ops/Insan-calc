import ctypes
import os
import sys


class C_CalcResult(ctypes.Structure):
  _fields_ = [
      ("total_cut_length_mm", ctypes.c_double),
      ("sheet_thickness_mm", ctypes.c_double),
      ("sheet_area_mm2", ctypes.c_double),
      ("pierce_count", ctypes.c_int),
      ("closed_loops", ctypes.c_int),
      ("open_chains", ctypes.c_int),
      ("primitives_count", ctypes.c_int),
      ("status_code", ctypes.c_int),
  ]


lib_path = os.path.abspath("./libinsan_calc.so")
if not os.path.exists(lib_path):
  print("Error: libinsan_calc.so not found. Run make first.")
  sys.exit(1)

insan_lib = ctypes.CDLL(lib_path)
insan_lib.insan_analyze_cad_file.argtypes = [ctypes.c_char_p]
insan_lib.insan_analyze_cad_file.restype = C_CalcResult


def ensure_test_dxf(path):
  if os.path.exists(path):
    return
  dxf_data = (
      "0\nSECTION\n2\nENTITIES\n"
      "0\nLINE\n10\n0.0\n20\n0.0\n11\n100.0\n21\n0.0\n"
      "0\nCIRCLE\n10\n50.0\n20\n50.0\n40\n50.0\n"
      "0\nENDSEC\n0\nEOF\n"
  )
  with open(path, "w", encoding="utf-8") as f:
    f.write(dxf_data)


def analyze_cad(filepath):
  res = insan_lib.insan_analyze_cad_file(filepath.encode("utf-8"))
  if res.status_code != 0:
    raise RuntimeError("C++ kernel failed to parse file.")
  return {
      "Cut Length (mm)": round(res.total_cut_length_mm, 4),
      "Pierces (B0)": res.pierce_count,
      "Closed Loops": res.closed_loops,
      "Open Chains": res.open_chains,
      "Primitives": res.primitives_count,
      "Thickness (mm)": round(res.sheet_thickness_mm, 2),
  }


if __name__ == "__main__":
  target = sys.argv[1] if len(sys.argv) > 1 else "test.dxf"
  if target == "test.dxf":
    ensure_test_dxf(target)

  stats = analyze_cad(target)
  print(f"=== INSAN C++/PYTHON BRIDGE: {target} ===")
  for k, v in stats.items():
    print(f"{k:18}: {v}")

  price = (stats["Cut Length (mm)"] / 1000.0) * 2.50 + stats[
      "Pierces (B0)"
  ] * 0.15
  print(f"Estimated Price   : ${price:.2f}")
