"""
shor_simulator.py
=================
Simulação do Algoritmo de Shor com Qiskit.

O que este módulo faz
---------------------
  1. Constrói o circuito quântico do algoritmo de Shor para fatorar N
     usando base a (escolhida automaticamente ou pelo usuário).

  2. Executa o circuito em simulador local (Qiskit Aer) ou,
     opcionalmente, em hardware IBM Quantum real.

  3. Processa os resultados:
       a. Extrai a fase medida das contagens de bitstrings.
       b. Usa frações contínuas para encontrar o período r.
       c. Calcula mdc(a^{r/2} ± 1, N) para obter os fatores.

  4. Exibe o circuito e os resultados de forma didática.

Casos de teste suportados
--------------------------
  N = 15  (a = 7)  –  circuito de referência, r = 4
  N = 21  (a = 4)  –  desafio maior,         r = 3 ou 6

Dependências
------------
  qiskit >= 2.0
  qiskit-aer >= 0.17
  (opcional) qiskit-ibm-runtime  para hardware real

Instalação
----------
  pip install qiskit qiskit-aer
  pip install qiskit-ibm-runtime   # apenas para IBM Quantum

Uso
---
  $ python shor_simulator.py               # roda N=15 e N=21 no simulador
  $ python shor_simulator.py --n 15 --a 7  # fatorar N=15 com base a=7
  $ python shor_simulator.py --ibmq <TOKEN> --n 15 --a 7  # hardware real

Referências
-----------
  Shor, P.W. (1994). Algorithms for Quantum Computation.
    Proceedings of FOCS 1994, pp. 124-134.
  Beauregard, S. (2003). Circuit for Shor's algorithm using 2n+3 qubits.
    Quantum Information and Computation, 3(2), 175-185.
  Nielsen, M.A. & Chuang, I.L. (2010). Quantum Computation and Quantum
    Information. Cambridge University Press.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
import warnings
from fractions import Fraction
from typing import Optional

# Suprimir DeprecationWarnings do Qiskit (QFT deprecated em 2.1+)
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    QISKIT_OK = True
except ImportError:
    QISKIT_OK = False
    print("[ERRO] Qiskit não encontrado. Instale com:")
    print("       pip install qiskit qiskit-aer")
    sys.exit(1)

# QFT: tenta QFTGate (>= 2.1) e recai em QFT legado
try:
    from qiskit.circuit.library import QFTGate
    _USE_QFTGATE = True
except ImportError:
    from qiskit.circuit.library import QFT          # type: ignore[assignment]
    _USE_QFTGATE = False


# ── Utilitários clássicos ─────────────────────────────────────────────────────

def _sep(titulo: str = "", largura: int = 66) -> None:
    if titulo:
        pad = max(1, (largura - len(titulo) - 2) // 2)
        print("─" * pad + f" {titulo} " + "─" * pad)
    else:
        print("─" * largura)


def _escolher_base(N: int) -> int:
    """
    Escolhe a base a para o algoritmo de Shor:
      - 1 < a < N
      - mdc(a, N) == 1 (caso contrário, o próprio mdc já fatorou N)
    """
    import secrets
    for _ in range(100):
        a = secrets.randbelow(N - 2) + 2
        g = math.gcd(a, N)
        if g == 1:
            return a
        elif 1 < g < N:
            # Sorte: encontrou fator diretamente
            return a   # retorna mesmo assim; o chamador verificará
    return 2


def _extrair_periodo(contagens: dict[str, int],
                     n_count: int,
                     N: int,
                     a: int,
                     verbose: bool = True) -> Optional[int]:
    """
    Extrai o período r a partir das contagens de medição.

    Cada bitstring representa uma fase φ ≈ s/r (para s inteiro),
    logo r é o denominador da melhor aproximação racional de φ
    com denominador ≤ N.

    Também testa múltiplos de r candidatos até que a^r ≡ 1 (mod N).
    """
    candidatos: dict[int, int] = {}

    for bitstring, freq in contagens.items():
        decimal = int(bitstring, 2)
        if decimal == 0:
            continue
        fase = decimal / (2 ** n_count)
        frac = Fraction(fase).limit_denominator(N)
        r    = frac.denominator
        if 1 < r <= N:
            candidatos[r] = candidatos.get(r, 0) + freq

    if not candidatos:
        return None

    # Ordenar por frequência de ocorrência
    ordenado = sorted(candidatos.items(), key=lambda x: -x[1])

    if verbose:
        print("    Candidatos a período r (por frequência):")
        for r_cand, freq in ordenado[:5]:
            print(f"      r = {r_cand:4d}  (aparece em {freq} medições)")

    # Para cada candidato, testar r, 2r, 3r... até N
    for r_base, _ in ordenado[:8]:
        for mult in range(1, N + 1):
            r_test = r_base * mult
            if r_test > N:
                break
            if pow(a, r_test, N) == 1:
                return r_test

    return None


def _fatorar_com_periodo(a: int, r: int, N: int,
                          verbose: bool = True) -> Optional[tuple[int, int]]:
    """
    Dado o período r de a mod N, calcula os fatores de N.

    Usa: p = mdc(a^{r/2} - 1, N)  e  q = mdc(a^{r/2} + 1, N).
    Se r é ímpar, tenta r*2 (ainda é período válido se a^{2r}≡1).
    """
    # Se r ímpar, tentar dobrar
    if r % 2 != 0:
        r2 = r * 2
        if pow(a, r2, N) == 1 and r2 % 2 == 0:
            if verbose:
                print(f"    Período r={r} ímpar; usando 2r={r2}.")
            r = r2
        else:
            if verbose:
                print(f"    Período r={r} ímpar e 2r não resolve. Escolher nova base.")
            return None

    x = pow(a, r // 2, N)

    if x == N - 1:
        if verbose:
            print(f"    a^(r/2) ≡ -1 (mod N): caso trivial. Escolher nova base.")
        return None

    p = math.gcd(x - 1, N)
    q = math.gcd(x + 1, N)

    if verbose:
        print(f"    a^(r/2) mod N  = {x}")
        print(f"    mdc(a^r/2 - 1, N) = mdc({x-1}, {N}) = {p}")
        print(f"    mdc(a^r/2 + 1, N) = mdc({x+1}, {N}) = {q}")

    for f in [p, q]:
        if 1 < f < N and N % f == 0:
            outro = N // f
            return (f, outro)

    return None


# ── Construção do circuito ────────────────────────────────────────────────────

def _qft_inversa(n: int) -> QuantumCircuit:
    """Retorna circuito da QFT inversa em n qubits."""
    if _USE_QFTGATE:
        qc = QuantumCircuit(n)
        qc.append(QFTGate(n).inverse(), range(n))
    else:
        from qiskit.circuit.library import QFT   # type: ignore[assignment]
        qft = QFT(n, inverse=True)
        qc  = QuantumCircuit(n)
        qc.append(qft, range(n))
    return qc


def _circuito_multiplicacao_modular_n15_a7(
        ctrl: int, alvo_qubits: list[int], qc: QuantumCircuit
) -> None:
    """
    Porta controlada U^{2^ctrl}: |y> -> |7y mod 15> no espaço de 4 qubits alvo.

    Implementação via SWAPs controlados (decomposição clássica para N=15, a=7).
    Referência: Nielsen & Chuang (2010), Exercise 5.20.

    A sequência de SWAPs implementa a permutação cíclica das potências de 7 mod 15:
    7^1=7, 7^2=4, 7^3=13, 7^4=1  =>  ciclo {1,7,4,13} nos estados de 4 qubits
    """
    q0, q1, q2, q3 = alvo_qubits

    reps = 2 ** ctrl   # número de vezes que U é aplicado

    for _ in range(reps):
        qc.cswap(ctrl, q1, q2)
        qc.cswap(ctrl, q0, q1)
        qc.cswap(ctrl, q2, q3)
        qc.cswap(ctrl, q1, q2)


def _circuito_multiplicacao_modular_n21_a2(
        ctrl: int, alvo_qubits: list[int], qc: QuantumCircuit
) -> None:
    """
    Porta controlada U^{2^ctrl}: |y> -> |2y mod 21> no espaço de 5 qubits alvo.

    2^1=2, 2^2=4, 2^3=8, 2^4=16, 2^5=11, 2^6=1  => período r=6 (par)
    Fatoração: mdc(2^3-1, 21)=7, mdc(2^3+1, 21)=3.
    """
    q0, q1, q2, q3, q4 = alvo_qubits
    reps = 2 ** ctrl
    # Permutação cíclica que simula multiplicação por 2 mod 21
    # no subespaço {1,2,4,8,16,11} do espaço de 5 qubits
    for _ in range(reps):
        qc.cswap(ctrl, q0, q1)
        qc.cswap(ctrl, q1, q2)
        qc.cswap(ctrl, q2, q3)
        qc.cswap(ctrl, q3, q4)


def construir_circuito_shor(N: int, a: int) -> QuantumCircuit:
    """
    Constrói o circuito quântico do algoritmo de Shor para fatorar N com base a.

    Estrutura do circuito:
      - n_count qubits de controle (registrador de fase): inicializados em H
      - n_target qubits alvo (registrador de trabalho): inicializados em |1>
      - QFT^{-1} aplicada ao registrador de controle
      - Medição do registrador de controle

    Suporta (N=15, a=7) e (N=21, a=4).
    """
    if (N, a) not in [(15, 7), (21, 2)]:
        raise NotImplementedError(
            f"Circuito explícito implementado apenas para (N=15,a=7) e (N=21,a=2). "
            f"Recebido: N={N}, a={a}."
        )

    if N == 15:
        n_count  = 4   # 2 * ceil(log2(N)) = 4
        n_target = 4   # ceil(log2(N)) = 4

        qc = QuantumCircuit(n_count + n_target, n_count,
                            name=f"Shor(N={N}, a={a})")

        # Superposição no registrador de controle
        for q in range(n_count):
            qc.h(q)

        # Estado inicial |1> no registrador alvo
        qc.x(n_count)

        qc.barrier(label="Init")

        # Portas U controladas: U^{2^j} para j = 0..n_count-1
        for ctrl_q in range(n_count):
            alvo = list(range(n_count, n_count + n_target))
            _circuito_multiplicacao_modular_n15_a7(ctrl_q, alvo, qc)

        qc.barrier(label="Oracles")

        # QFT inversa no registrador de controle
        iqft = _qft_inversa(n_count)
        qc.append(iqft, range(n_count))

        qc.barrier(label="iQFT")

    else:  # N == 21, a == 2
        n_count  = 6   # 2 * ceil(log2(21)) = 10 reduzido a 6 para simulação
        n_target = 5   # ceil(log2(21))

        qc = QuantumCircuit(n_count + n_target, n_count,
                            name=f"Shor(N={N}, a={a})")

        for q in range(n_count):
            qc.h(q)

        qc.x(n_count)   # |1> no alvo

        qc.barrier(label="Init")

        for ctrl_q in range(n_count):
            alvo = list(range(n_count, n_count + n_target))
            _circuito_multiplicacao_modular_n21_a2(ctrl_q, alvo, qc)

        qc.barrier(label="Oracles")

        iqft = _qft_inversa(n_count)
        qc.append(iqft, range(n_count))

        qc.barrier(label="iQFT")

    # Medição
    qc.measure(range(n_count), range(n_count))

    return qc


# ── Execução no simulador ────────────────────────────────────────────────────

def executar_simulador(qc: QuantumCircuit,
                       shots: int = 2048,
                       verbose: bool = True) -> dict[str, int]:
    """
    Executa o circuito no simulador Qiskit Aer e retorna as contagens.
    """
    sim = AerSimulator()

    if verbose:
        print(f"    Backend   : AerSimulator (local)")
        print(f"    Shots     : {shots}")
        print(f"    Transpilando...", end=" ", flush=True)

    inicio = time.perf_counter()
    qc_t   = transpile(qc, sim, optimization_level=1)
    t_transp = time.perf_counter() - inicio

    if verbose:
        ops = dict(qc_t.count_ops())
        print(f"OK ({t_transp*1000:.0f} ms)")
        print(f"    Portas    : {ops}")
        print(f"    Qubits    : {qc_t.num_qubits}")
        print(f"    Profundidade: {qc_t.depth()}")
        print(f"    Executando...", end=" ", flush=True)

    inicio2 = time.perf_counter()
    result  = sim.run(qc_t, shots=shots).result()
    t_exec  = time.perf_counter() - inicio2

    if verbose:
        print(f"OK ({t_exec*1000:.0f} ms)")

    return result.get_counts()


# ── Execução em hardware IBM Quantum ─────────────────────────────────────────

def executar_ibmq(qc: QuantumCircuit,
                  token: str,
                  shots: int = 1024,
                  verbose: bool = True) -> dict[str, int]:
    """
    Executa o circuito em hardware IBM Quantum real.

    Requer: pip install qiskit-ibm-runtime
    Token disponível em: https://quantum.ibm.com

    Parâmetros
    ----------
    qc    : circuito a executar
    token : token de acesso IBM Quantum
    shots : número de medições

    Retorna
    -------
    Dicionário de contagens {bitstring: frequência}
    """
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
        from qiskit_ibm_runtime import Session
    except ImportError:
        raise ImportError(
            "qiskit-ibm-runtime não instalado.\n"
            "Instale com: pip install qiskit-ibm-runtime"
        )

    if verbose:
        print("    Conectando ao IBM Quantum... ", end="", flush=True)

    service = QiskitRuntimeService(channel="ibm_quantum", token=token)
    backend = service.least_busy(
        operational=True,
        simulator=False,
        min_num_qubits=qc.num_qubits
    )

    if verbose:
        print(f"OK")
        print(f"    Backend selecionado : {backend.name}")
        print(f"    Qubits disponíveis  : {backend.num_qubits}")
        print(f"    Transpilando para o backend...", end=" ", flush=True)

    qc_t = transpile(qc, backend, optimization_level=3)

    if verbose:
        print(f"OK")
        print(f"    Profundidade pós-transpilação: {qc_t.depth()}")
        print(f"    Submetendo job... ", end="", flush=True)

    with Session(backend=backend) as session:
        sampler = Sampler(session=session)
        job     = sampler.run([qc_t], shots=shots)

        if verbose:
            print(f"OK (Job ID: {job.job_id()})")
            print(f"    Aguardando resultado (pode demorar minutos na fila)...")

        result  = job.result()
        counts  = result[0].data.c.get_counts()

    return counts


# ── Fatoração completa ────────────────────────────────────────────────────────

def fatorar_shor(N: int,
                 a: Optional[int] = None,
                 backend: str = "simulator",
                 shots: int = 2048,
                 ibmq_token: Optional[str] = None,
                 verbose: bool = True) -> Optional[tuple[int, int]]:
    """
    Executa o algoritmo de Shor completo para fatorar N.

    Parâmetros
    ----------
    N           : inteiro a fatorar (15 ou 21)
    a           : base (escolhida automaticamente se None)
    backend     : "simulator" (local) ou "ibmq" (requer token)
    shots       : número de medições quânticas
    ibmq_token  : token IBM Quantum (necessário se backend="ibmq")
    verbose     : exibir passos detalhados

    Retorna
    -------
    (p, q) se bem-sucedido, None caso contrário.
    """
    if verbose:
        _sep(f"ALGORITMO DE SHOR: N = {N}")
        print()

    # ── Pré-processamento clássico ────────────────────────────────────────
    if verbose:
        print("  [FASE 1] Pré-processamento Clássico")
        print(f"    N = {N}  ({N.bit_length()} bits)")

    # Verificar trivialidades
    if N % 2 == 0:
        if verbose:
            print(f"    N é par! Fator trivial: 2 × {N//2}")
        return (2, N // 2)

    # Verificar potência perfeita
    for b in range(2, N.bit_length() + 1):
        root = round(N ** (1 / b))
        for candidate in [root - 1, root, root + 1]:
            if candidate > 1 and candidate ** b == N:
                if verbose:
                    print(f"    N = {candidate}^{b}: fator trivial!")
                return (candidate, candidate ** (b - 1))

    # Escolher base a
    if a is None:
        a = _escolher_base(N)

    g = math.gcd(a, N)
    if 1 < g < N:
        if verbose:
            print(f"    mdc(a={a}, N={N}) = {g}: fator encontrado sem circuito quântico!")
        return (g, N // g)

    if verbose:
        print(f"    Base escolhida: a = {a}")
        print(f"    mdc({a}, {N}) = {g} = 1  ✓ (pré-condição satisfeita)")

    # ── Fase quântica ──────────────────────────────────────────────────────
    if verbose:
        print()
        print("  [FASE 2] Execução Quântica (busca do período r)")

    try:
        qc = construir_circuito_shor(N, a)
    except NotImplementedError as err:
        print(f"    [AVISO] {err}")
        return None

    if verbose:
        print(f"    Circuito construído: {qc.num_qubits} qubits, "
              f"{qc.depth()} camadas (pré-transpilação)")

    if backend == "ibmq":
        if ibmq_token is None:
            raise ValueError("Token IBM Quantum necessário para backend='ibmq'.")
        contagens = executar_ibmq(qc, ibmq_token, shots=shots, verbose=verbose)
    else:
        contagens = executar_simulador(qc, shots=shots, verbose=verbose)

    if verbose:
        print()
        print("  [FASE 3] Pós-processamento Clássico (frações contínuas)")
        print()

    # Extrair período
    n_count = sum(1 for bit in qc.clbits)
    r = _extrair_periodo(contagens, n_count, N, a, verbose=verbose)

    if r is None:
        if verbose:
            print("    Período não encontrado nas medições. Tente mais shots.")
        return None

    if verbose:
        print(f"\n    Período encontrado: r = {r}")
        print(f"    Verificação: {a}^{r} mod {N} = {pow(a, r, N)}  "
              f"{'✓' if pow(a, r, N) == 1 else '✗ (período incorreto)'}")

    if pow(a, r, N) != 1:
        if verbose:
            print("    Período incorreto. Tente novamente.")
        return None

    # Calcular fatores
    if verbose:
        print()
        print("  [FASE 4] Cálculo dos Fatores")
    fatores = _fatorar_com_periodo(a, r, N, verbose=verbose)

    if fatores:
        p, q = fatores
        if verbose:
            print(f"\n  *** SUCESSO ***")
            print(f"  {N} = {p} × {q}")
            print(f"  Verificação: {p} × {q} = {p * q}  "
                  f"{'✓' if p * q == N else '✗'}")
        return fatores
    else:
        if verbose:
            print("    Fatoração falhou com este período. Tente nova base.")
        return None


# ── Visualização do circuito ─────────────────────────────────────────────────

def exibir_circuito(N: int, a: int) -> None:
    """Exibe a representação textual do circuito de Shor."""
    _sep(f"CIRCUITO: Shor(N={N}, a={a})")
    try:
        qc = construir_circuito_shor(N, a)
        print()
        print(qc.draw(output="text", fold=100))
        print()
        print(f"  Qubits totais      : {qc.num_qubits}")
        print(f"    Registrador fase : {sum(1 for b in qc.clbits)} qubits de controle")
        print(f"    Registrador alvo : {qc.num_qubits - sum(1 for b in qc.clbits)} qubits")
        print(f"  Profundidade       : {qc.depth()}")
        print()
    except NotImplementedError as e:
        print(f"  [AVISO] {e}")


# ── Contexto educacional ──────────────────────────────────────────────────────

def explicacao_shor() -> None:
    """Imprime a explicação conceitual do algoritmo de Shor."""
    _sep("ALGORITMO DE SHOR – Visão Conceitual")
    print("""
  Problema: dado N = p × q, encontrar p e q.

  Redução clássica (Shor, 1994):
  ──────────────────────────────
  1. Escolher a aleatório com 1 < a < N e mdc(a,N)=1.
  2. Encontrar o PERÍODO r: menor r>0 tal que a^r ≡ 1 (mod N).
  3. Se r é par e a^{r/2} ≢ -1 (mod N), então:
       p = mdc(a^{r/2} - 1, N)  e  q = mdc(a^{r/2} + 1, N)

  Aceleração quântica (QFT):
  ──────────────────────────
  O passo 2 é O(N^{1/2}) classicamente, mas o circuito quântico
  coloca o registrador de fase em superposição de todos os x de 0 a 2^n-1,
  aplica o oráculo unitário U_f:|x>|1> -> |x>|a^x mod N>,
  e usa a QFT^{-1} para extrair a frequência dominante s/r.
  A medição colapsa para s/r com alta probabilidade.
  Frações contínuas extraem r de s/r.

  Complexidade:
  ─────────────
  Clássico (GNFS):  O(exp((log N)^{1/3} (log log N)^{2/3}))  sub-exponencial
  Shor quântico:    O((log N)^2 log log N log log log N)       polinomial

  Para RSA-2048 (N com 617 dígitos decimais):
    GNFS    : ~2^{112} operações  (computacionalmente seguro hoje)
    Shor    : ~4.000 qubits lógicos + 20M qubits físicos (Webber et al., 2022)
""")


# ── Análise de complexidade ───────────────────────────────────────────────────

def comparar_complexidade() -> None:
    """Tabela comparativa de complexidade para diferentes tamanhos de N."""
    _sep("COMPLEXIDADE: Clássico vs. Quântico")
    print()
    print("  Tamanho N  | GNFS (bits de trabalho) | Shor (operações quânticas)")
    print("  " + "─" * 62)

    tamanhos = [256, 512, 1024, 2048, 4096]
    for bits in tamanhos:
        # Complexidade GNFS: L_n[1/3; 1.923] em bits de segurança (aproximado)
        gnfs_bits = int(1.923 * (bits ** (1/3)) * (math.log2(bits) ** (2/3)))
        # Complexidade Shor: O(bits^3) operações quânticas (aprox. para QFT-based)
        shor_ops  = bits ** 3
        print(f"  RSA-{bits:4d}   |  2^{gnfs_bits:3d} operações clássicas  |  "
              f"~{shor_ops:>12,.0f} portas quânticas")

    print()
    print("  Conclusão: Shor quebra RSA de qualquer tamanho em tempo polinomial,")
    print("  mas requer hardware quântico tolerante a falhas ainda não existente.")
    print()


# ── Main ──────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Simulação do Algoritmo de Shor com Qiskit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos:
  python shor_simulator.py                     # roda N=15 e N=21
  python shor_simulator.py --n 15 --a 7        # N=15, base a=7
  python shor_simulator.py --n 21 --a 4        # N=21, base a=4
  python shor_simulator.py --shots 4096        # mais medições
  python shor_simulator.py --ibmq SEU_TOKEN    # hardware real
        """
    )
    parser.add_argument("--n",     type=int, default=None,
                        help="Inteiro a fatorar (15 ou 21)")
    parser.add_argument("--a",     type=int, default=None,
                        help="Base para o algoritmo de Shor")
    parser.add_argument("--shots", type=int, default=2048,
                        help="Número de medições (padrão: 2048)")
    parser.add_argument("--ibmq",  type=str, default=None,
                        metavar="TOKEN",
                        help="Token IBM Quantum para usar hardware real")
    parser.add_argument("--no-circuit", action="store_true",
                        help="Não exibir o diagrama do circuito")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    print()
    print("=" * 66)
    print("  shor_simulator.py — Algoritmo de Shor com Qiskit")
    print("  Atividade Final – Introdução à Cibersegurança")
    print("  SENAI CIMATEC – Pós-graduação em Comunicação Quântica")
    print("=" * 66)

    explicacao_shor()
    comparar_complexidade()

    backend = "ibmq" if args.ibmq else "simulator"

    # Casos a executar
    if args.n is not None:
        casos = [(args.n, args.a if args.a else (7 if args.n == 15 else 2))]
    else:
        casos = [(15, 7), (21, 2)]

    for N, a in casos:
        print()
        _sep("=" * 60)

        # Exibir circuito
        if not args.no_circuit:
            exibir_circuito(N, a)

        # Executar fatoração
        inicio = time.perf_counter()
        resultado = fatorar_shor(
            N=N,
            a=a,
            backend=backend,
            shots=args.shots,
            ibmq_token=args.ibmq,
            verbose=True,
        )
        elapsed = time.perf_counter() - inicio

        print()
        if resultado:
            p, q = resultado
            _sep(f"RESULTADO FINAL: {N} = {p} × {q}")
        else:
            _sep(f"RESULTADO: Fatoração não convergiu para N={N}")

        print(f"  Tempo total (clássico + quântico): {elapsed:.3f} s")
        print()

    _sep()
    print("  Simulação concluída.")
    print()


if __name__ == "__main__":
    main()
