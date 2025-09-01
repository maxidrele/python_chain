import argparse, json, os
import numpy as np
from qutip import tensor, qeye, destroy, mesolve, Options

# Función auxiliar para crear los operadores (la que tenías antes)
def cadena(dim_resonador: int, n_transmons: int, dim_transmon: int, return_transmons: bool = True):
    a_r1 = tensor(destroy(dim_resonador),
                  *(qeye(dim_transmon) for _ in range(n_transmons)),
                  qeye(dim_resonador))

    as_t = []
    if return_transmons:

      for n in range(n_transmons):
          at_n = tensor(
              qeye(dim_resonador),
              *(qeye(dim_transmon) for _ in range(n)),
              destroy(dim_transmon),
              *(qeye(dim_transmon) for _ in range(n + 1, n_transmons)),
              qeye(dim_resonador)
          )
          as_t.append(at_n)

      a_rN = tensor(qeye(dim_resonador),
                    *(qeye(dim_transmon) for _ in range(n_transmons)),
                    destroy(dim_resonador))

    return [a_r1] + as_t + [a_rN]


def hamiltoniano(ops, params: dict):
    """
    Args:
        n_transmons: Número de transmones en la cadena
        params: Diccionario con los parámetros del sistema que debe contener:
            - omega_r1: Frecuencia resonador izquierdo
            - omega_rN: Frecuencia resonador derecho
            - omega_t: Lista de frecuencias de los transmones [ω_t1, ω_t2,...]
            - alpha_t: Lista de anarmonicidades de los transmones [α_t1, α_t2,...]
            - g_t: Lista de acoplamientos entre transmones vecinos [g_t1_t2, g_t2_t3,...]
            - g_r1_t1: Acoplamiento resonador izquierdo-primer transmon
            - g_tN_rN: Acoplamiento último transmon-resonador derecho

    Returns:
        El Hamiltoniano total del sistema
    """
    # Primero creamos los operadores de destrucción
    a_r1, *a_ts, a_rN = ops
    n_transmons = len(a_ts)

    # Operadores de creación
    a_r1_dag = a_r1.dag()
    a_rN_dag = a_rN.dag()
    a_ts_dag = [a.dag() for a in a_ts]

    # Inicializamos el Hamiltoniano
    H = 0

    # Términos de los resonadores
    H += params['w_r1'] * a_r1_dag * a_r1

    # Acoplamiento entre el primer resonador y su vecino
    H += +params['g_r1_t1'] * (a_r1 * a_ts_dag[0] + a_r1_dag * a_ts[0])

    # Términos de los transmones
    for i in range(n_transmons):
        a_t = a_ts[i]
        a_t_dag = a_ts_dag[i]

        # Términos de energía y anarmonicidad
        H += params['w_t'] * a_t_dag * a_t
        H += -params['alpha_t'] * a_t_dag * a_t_dag * a_t * a_t / 2

        # Acoplamientos entre transmones vecinos
        if i < n_transmons - 1:
            H += params['g_t'] * (a_t * a_ts_dag[i+1] + a_t_dag * a_ts[i+1])

    # Acoplamiento del ultimo transmon con el resonador
    H += +params['g_tN_rN'] * (a_ts[-1] * a_rN_dag + a_ts_dag[-1] * a_rN)

    # Hamiltoniano del ultimo resonador
    H += params['w_rN'] * a_rN_dag * a_rN
    return H


def jumps(a_r1, a_rN, gamma_1, gamma_N, n_r1, n_rN):
    b_r1 = np.sqrt(gamma_1*(n_r1+1))*a_r1
    b_r1_dag = np.sqrt(gamma_1*n_r1)*a_r1.dag()
    b_rN = np.sqrt(gamma_N*(n_rN+1))*a_rN
    b_rN_dag = np.sqrt(gamma_N*n_rN)*a_rN.dag()

    return [b_r1, b_r1_dag, b_rN, b_rN_dag]

def simular(dim_resonador, n_transmons, dim_transmon, H_params, jumps_params, sim_params, psi_0, return_obs = True):

  a_s = cadena(dim_resonador, n_transmons, dim_transmon) #creamos los ops de destrucción de la cadena

  H = hamiltoniano(a_s, H_params) #con esto definimos el Hamiltoniano

  jumps_ops = jumps(a_s[0], a_s[-1], gamma_1=jumps_params["gamma_1"], gamma_N=jumps_params["gamma_N"],
                    n_r1=jumps_params["n_h"], n_rN=jumps_params["n_c"])

  times = np.linspace(0.0, sim_params["T_final"], sim_params["nsteps"])

  options = {'progress_bar': True}

  e_ops_list = []
  if return_obs:
    e_ops_list = [a.dag()*a for a in a_s]
    e_ops_list.append((a_s[0].dag()*a_s[1] + a_s[1].dag()*a_s[0]))
 #correlacion entre el primer reson


  result = mesolve(H, psi_0, times, jumps_ops, e_ops=e_ops_list, options=options)
  #result = mesolve(H, psi_0, times, jumps_ops, e_ops=e_ops_list)

  return times, result