#!/usr/bin/env python3
"""Evaluate a GI-disk reference state with Athena's compiled Helmholtz EOS."""

import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / (
    "inputs/hydro/"
    "athinput.gi_disk_v42_model5_q0p5_ye0p5_hatc0p05_"
    "rout200_r240_diskmonopole_diode_thetacut0p3_1omega_4rank"
)
C = 2.99792458e10
G = 6.67430e-8
MSUN = 1.98847e33
KB = 1.38e-16
HBAR = 1.054e-27
MP = 1.67e-24
NA = 6.02214076e23
MEV_TO_ERG = 1.602176634e-6


def fermi3(eta):
    a = math.exp(-abs(eta))
    s = eta**4 / 4.0 + 4.9348022 * eta**2 + 11.351273
    ff = 6.0 * (a - a**2 / 16.0 + a**3 / 81.0 - a**4 / 256.0)
    return ff if eta < 0.0 else s - ff


def fermi4(eta):
    a = math.exp(-abs(eta))
    s = eta**5 / 5.0 + 6.5797363 * eta**3 + 45.457576 * eta
    ff = 24.0 * (a - a**2 / 32.0 + a**3 / 243.0)
    return ff if eta < 0.0 else s + ff


def fermi5(eta):
    a = math.exp(-abs(eta))
    s = eta**6 / 6.0 + 8.2246703 * eta**4 + 113.64394 * eta**2 + 236.53226
    ff = 120.0 * (a - a**2 / 64.0 + a**3 / 729.0)
    return ff if eta < 0.0 else s - ff


def pgen_eta(rho, temp, ye):
    third = 1.0 / 3.0
    a = (KB / (HBAR * C))**3 * temp**3 * (MP / rho) * (math.pi / 3.0)
    b = math.sqrt(4.0 * a**6 + 27.0 * a**4 * ye**2)
    term = (9.0 * a**2 * ye + math.sqrt(3.0) * b)**third
    norm = 6.0**(2.0 * third)
    return math.pi * (2.0**third * term / (norm * a)
                      - 2.0 * 3.0**third * a / (norm * term))


def cooling(rho, temp, ye, xalpha):
    eta = pgen_eta(rho, temp, ye)
    temp_mev = temp * 8.617333262e-11
    yp = max(0.0, ye - 0.5 * xalpha)
    yn = max(0.0, 1.0 - ye - 0.5 * xalpha)
    capture = (2.073 * NA * temp_mev**6 * MEV_TO_ERG
               * (yp * fermi5(eta) + yn * fermi5(-eta)) / 118.266)
    pair_factor = (fermi4(eta) * fermi3(-eta)
                   + fermi4(-eta) * fermi3(eta)) / (2.0 * 23.3309 * 5.6822)
    pair = (0.145 * NA * 1.0e8 * temp_mev**9 / rho
            * pair_factor * MEV_TO_ERG)
    return eta, yp, yn, capture, pair


def eos_query(rho, pressure, ye, input_path, ranks):
    text = input_path.read_text()
    probe = ("<problem>\n"
             "eos_reference_probe = true\n"
             f"eos_probe_rho = {rho:.17e}\n"
             f"eos_probe_pressure = {pressure:.17e}\n"
             f"eos_probe_ye = {ye:.17e}\n")
    if "<problem>" not in text:
        raise RuntimeError(f"No <problem> block in {input_path}")
    text = text.replace("<problem>\n", probe, 1)
    with tempfile.TemporaryDirectory(prefix="athena-reference-") as tmp:
        tmpdir = Path(tmp)
        generated = tmpdir / "athinput.reference"
        generated.write_text(text)
        athena = [str(ROOT / "bin/athena"), "-i", str(generated),
                  f"job/problem_id={tmpdir / 'probe'}", "time/nlim=0"]
        if ranks == 1:
            # Spherical-harmonic gravity requires one MeshBlock per rank.
            athena += ["meshblock/nx1=20", "meshblock/nx2=48", "meshblock/nx3=16"]
        command = athena if ranks == 1 else ["mpirun", "-np", str(ranks)] + athena
        run = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        output = run.stdout + run.stderr
        if run.returncode != 0:
            raise RuntimeError(f"Athena EOS probe failed:\n{output}")
    match = re.search(r"^EOS_REFERENCE\s+(.+)$", output, re.MULTILINE)
    if not match:
        raise RuntimeError(
            "Athena output did not contain EOS_REFERENCE; rebuild bin/athena.\n"
            + output
        )
    return {key: float(value) for key, value in re.findall(r"(\w+)=([^\s]+)", match.group(1))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mass-msun", type=float, required=True)
    parser.add_argument("--ye", type=float, required=True)
    parser.add_argument("--chat", type=float, required=True,
                        help="initialization proxy sqrt(P/rho)/c")
    parser.add_argument("--q", type=float, required=True, help="initial Toomre Q")
    parser.add_argument("--radius-rg", type=float, default=100.0)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    mass = args.mass_msun * MSUN
    omega = C**3 / (G * mass * args.radius_rg**1.5)
    rho = omega**2 / (math.pi * G * args.q * math.sqrt(2.0 * math.pi))
    pressure = rho * (args.chat * C)**2
    eos = eos_query(rho, pressure, args.ye, args.input.resolve(), args.ranks)
    eta, yp, yn, capture, pair = cooling(
        rho, eos["temperature"], args.ye, eos["xalpha"])
    total = capture + pair
    tcool = eos["energy_density"] / (rho * total)
    result = {
        "mass_msun": args.mass_msun, "radius_rg": args.radius_rg,
        "q_initial": args.q, "ye": args.ye, "chat_initial": args.chat,
        "omega_s-1": omega, "omega_inverse_s": 1.0 / omega,
        "rho_g_cm3": rho, "pressure_erg_cm3": pressure,
        "temperature_K": eos["temperature"],
        "energy_density_erg_cm3": eos["energy_density"],
        "cs_over_c_eos": eos["cs_over_c"], "xalpha": eos["xalpha"],
        "eta": eta, "yp_free": yp, "yn_free": yn,
        "capture_cooling_erg_g_s": capture,
        "pair_cooling_erg_g_s": pair,
        "total_cooling_erg_g_s": total,
        "pair_fraction": pair / total, "cooling_time_s": tcool,
        "omega_tcool": omega * tcool,
    }
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        for key, value in result.items():
            print(f"{key:32s} {value:.9e}" if isinstance(value, float) else f"{key:32s} {value}")


if __name__ == "__main__":
    main()
