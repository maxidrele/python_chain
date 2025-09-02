import numpy as np
from qutip import tensor, basis
# Importar las funciones del otro archivo
from simular_cadena import cadena, hamiltoniano, jumps, simular

n_transmons = 3

omega = 5

omega_h = 3  # Temperatura de las fuentes (2-8 GHz = 380 - 96 K)
omega_c = 2

g = 0.02 # 0.02 GHz (20 MHz)

alpha_t = 0.1

H_params = {
    'w_r1': omega,
    'w_rN': omega,
    'w_t': omega,  # frecuencias de los transmones
    'alpha_t': alpha_t,  # anarmonicidades
    'g_t': g,          # acoplamiento t_i-t_i+1
    'g_r1_t1': g,        # acoplamiento r1-t1
    'g_tN_rN': g        # acoplamiento tN-rN
  }

jumps_params = {
    'gamma_1': 0.5/(2*np.pi),     #gamma_1
    'gamma_N': 0.5/(2*np.pi),     #gamma_N
    'n_h': 1/(np.exp(omega/omega_h)-1),
    'n_c': 1/(np.exp(omega/omega_c)-1)
  }

sim_params = {
    'T_final': 10000,
    'nsteps' : 1000
  }

d=4  # Dimensión del espacio de Hilbert de cada elemento (resonadores y transmones)
psi_0 = tensor(basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 0), basis(d, 0))  # Estado inicial: todos en el estado base   .

times, result = simular(dim_resonador=d, n_transmons=n_transmons, dim_transmon=d, H_params=H_params, jumps_params=jumps_params, sim_params=sim_params, psi_0= psi_0)

import json
from datetime import datetime

print("✓ Simulación completada!")
# 1. Guardar en formato NPZ (recomendado para numpy arrays)
np.savez('resultados_simulacion_d4.npz', 
         times=times,
         expect=result.expect,
         params=H_params,  # Guardar también los parámetros
         jumps_params=jumps_params,
         sim_params=sim_params)

print("✓ Datos guardados en 'resultados_simulacion_d4.npz'")

# 2. Guardar en formato CSV (para Excel o otros programas)


# 3. Guardar metadatos en JSON (para referencia)
metadata = {
    'fecha_simulacion': datetime.now().isoformat(),
    'parametros_hamiltoniano': H_params,
    'parametros_salto': jumps_params,
    'parametros_simulacion': sim_params,
    'n_transmons': n_transmons,
    'dimension': d,
    'numero_observables': len(result.expect),
    'descripcion_observables': [
        'Resonador 1: a†a',
        'Transmon 1: a†a',
        'Transmon 2: a†a', 
        'Transmon 3: a†a',
        'Resonador N: a†a',
        'Correlación resonador-transmon'
    ]
}

with open('metadata.json_d4', 'w') as f:
    json.dump(metadata, f, indent=4)

print("✓ Metadatos guardados en 'metadata.json'")

# 4. Crear un archivo de resumen
with open('resumen_simulacion_d4.txt', 'w') as f:
    f.write(f"RESUMEN DE SIMULACIÓN - {datetime.now()}\n")
    f.write("=" * 50 + "\n\n")
    f.write(f"Tiempo final: {sim_params['T_final']} unidades\n")
    f.write(f"Número de pasos: {sim_params['nsteps']}\n")
    f.write(f"Número de transmones: {n_transmons}\n")
    f.write(f"Dimensión: {d}\n")
    f.write(f"Número de observables: {len(result.expect)}\n\n")
    
    f.write("Valores finales de los observables:\n")
    for i, obs in enumerate(result.expect):
        f.write(f"  Observable {i+1}: {obs[-1]:.6f}\n")

print("✓ Resumen guardado en 'resumen_simulacion.txt'")