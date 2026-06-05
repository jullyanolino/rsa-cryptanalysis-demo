"""
rsa_attacks.py
==============
Simulação educacional dos principais ataques clássicos ao RSA.

Ataques implementados
---------------------
  1. Fatoração por Força Bruta
     - Testa todos os primos até sqrt(n)
     - Viável apenas para módulos pequenos (até ~30 bits)
     - Demonstra por que n deve ser grande

  2. Ataque de Wiener (expoente privado pequeno)
     - Explora frações contínuas quando d < n^{1/4}/3
     - Recupera d em tempo polinomial sem fatorar n
     - Referência: Wiener (1990)

  3. Ataque de Módulo Comum
     - Quando duas entidades compartilham n com expoentes diferentes
     - Recupera M sem fatorar n nem conhecer d
     - Demonstra por que módulos distintos são obrigatórios

  4. Ataque de Hastad (expoente público pequeno, múltiplos destinatários)
     - Quando a mesma mensagem M é enviada a e=3 destinatários
     - Usa o Teorema Chinês do Resto para recuperar M^3 em Z e depois raiz cúbica
     - Demonstra por que e pequeno + mesma mensagem sem padding é perigoso

Dependências
------------
  rsa_core.py (mesmo diretório)
  Biblioteca padrão Python >= 3.10 (math, itertools, secrets)

Uso
---
  $ python rsa_attacks.py

Referências
-----------
  Paar & Pelzl (2010), cap. 7 e 8
  Menezes, van Oorschot & Vanstone (2018), cap. 8
  Wiener, M. (1990). Cryptanalysis of Short RSA Secret Exponents.
    IEEE Transactions on Information Theory, 36(3), 553-558.
"""

from __future__ import annotations

import math
import secrets
import time
from typing import Optional

from rsa_core import (
    ChavePrivada, ChavePublica, ParDeChaves,
    cifrar_textbook, decifrar_textbook,
    euclides_estendido, exponenciacao_modular,
    gerar_primo, inverso_modular, mdc, miller_rabin,
)

# ── Utilitários ───────────────────────────────────────────────────────────────

def _sep(titulo: str = "", largura: int = 66) -> None:
    if titulo:
        pad = (largura - len(titulo) - 2) // 2
        print("─" * pad + f" {titulo} " + "─" * pad)
    else:
        print("─" * largura)


def _isqrt(n: int) -> int:
    """Raiz inteira por Newton-Raphson (sem math.isqrt para compatibilidade)."""
    if n < 0:
        raise ValueError("Raiz de número negativo.")
    if n == 0:
        return 0
    x = n
    y = (x + 1) // 2
    while y < x:
        x, y = y, (y + n // y) // 2
    return x


def _raiz_cubica_inteira(n: int) -> int:
    """Raiz cúbica inteira de n (Newton-Raphson)."""
    if n < 0:
        return -_raiz_cubica_inteira(-n)
    if n == 0:
        return 0
    x = int(round(n ** (1 / 3)))
    # Ajuste fino
    for delta in [-2, -1, 0, 1, 2]:
        c = x + delta
        if c >= 0 and c ** 3 == n:
            return c
    # Busca mais ampla se necessário
    lo, hi = 0, min(n, 2 ** ((n.bit_length() // 3) + 2))
    while lo <= hi:
        mid = (lo + hi) // 2
        c3  = mid ** 3
        if c3 == n:
            return mid
        elif c3 < n:
            lo = mid + 1
        else:
            hi = mid - 1
    return lo - 1


# ── Ataque 1: Fatoração por Força Bruta ──────────────────────────────────────

def fatorar_forca_bruta(n: int,
                         verbose: bool = True) -> Optional[tuple[int, int]]:
    """
    Fatora n por tentativa de divisão até sqrt(n).

    Retorna (p, q) com p <= q se bem-sucedido, None se n for primo.

    Complexidade: O(sqrt(n)) — intratável para n > 2^60 aprox.
    """
    if n < 4:
        return None

    limite = _isqrt(n) + 1
    candidatos = [2] + list(range(3, limite, 2))   # 2 e ímpares

    inicio = time.perf_counter()
    for p in candidatos:
        if n % p == 0:
            q = n // p
            elapsed = time.perf_counter() - inicio
            if verbose:
                print(f"    Tentativas testadas : até {p}")
                print(f"    Tempo decorrido     : {elapsed*1000:.3f} ms")
                print(f"    Fator encontrado    : p = {p}")
                print(f"    Cofator             : q = {q}")
                print(f"    Verificação         : p × q = {p * q} = n  ✓")
            return (p, q)

    elapsed = time.perf_counter() - inicio
    if verbose:
        print(f"    n parece primo (testados até sqrt(n)={limite}, em {elapsed*1000:.2f} ms)")
    return None


def demo_forca_bruta() -> None:
    """Demonstra a fatoração por força bruta em módulos RSA pequenos."""
    _sep("ATAQUE 1: Fatoração por Força Bruta")
    print("""
  Princípio: testar todos os inteiros até sqrt(n) como possíveis fatores.
  Funciona para n pequeno; torna-se astronomicamente lento para n >= 1024 bits.
""")

    casos = [
        ("Exemplo do relatório (p=61, q=53)", 61 * 53),
        ("Primos de 15 bits cada",  gerar_primo(15) * gerar_primo(15)),
        ("Primos de 18 bits cada",  gerar_primo(18) * gerar_primo(18)),
    ]

    for descricao, n in casos:
        print(f"  {descricao}  →  n = {n}  ({n.bit_length()} bits)")
        resultado = fatorar_forca_bruta(n, verbose=True)
        print()

    print("  Conclusão: a força bruta é inviável para n >= 64 bits.")
    print("  Para RSA-2048, testando 10^9 candidatos/s, levaria ~10^{290} anos.\n")


# ── Ataque 2: Wiener (expoente privado pequeno) ───────────────────────────────

def _convergentes(a: list[int]) -> list[tuple[int, int]]:
    """
    Calcula os convergentes da fração contínua [a0; a1, a2, ...].

    Retorna lista de pares (numerador, denominador).
    """
    convergentes = []
    h_prev, h_curr = 1, a[0]
    k_prev, k_curr = 0, 1

    convergentes.append((h_curr, k_curr))

    for i in range(1, len(a)):
        h_prev, h_curr = h_curr, a[i] * h_curr + h_prev
        k_prev, k_curr = k_curr, a[i] * k_curr + k_prev
        convergentes.append((h_curr, k_curr))

    return convergentes


def _fracao_continua(num: int, den: int) -> list[int]:
    """Expansão em fração contínua de num/den."""
    coeficientes = []
    while den:
        q = num // den
        coeficientes.append(q)
        num, den = den, num - q * den
    return coeficientes


def wiener_attack(e: int, n: int,
                  verbose: bool = True) -> Optional[int]:
    """
    Ataque de Wiener: recupera d quando d < n^{1/4} / 3.

    Baseia-se no fato de que k/d (onde e*d - k*phi(n) = 1)
    é um convergente da fração contínua de e/n.

    Retorna d se encontrado, None caso contrário.

    Referência: Wiener (1990), IEEE Transactions on Information Theory.
    """
    if verbose:
        print(f"    e  = {e}")
        print(f"    n  = {n}")
        print(f"    n^(1/4)/3 ≈ {int(n**0.25)//3}  (limiar de segurança)")

    coefs = _fracao_continua(e, n)
    convergentes = _convergentes(coefs)

    for k, d_cand in convergentes:
        if k == 0:
            continue
        # Se d_cand é o expoente privado real:
        # e * d_cand ≡ 1 (mod phi(n))  =>  e * d_cand = 1 + k * phi(n)
        # phi(n) = (e * d_cand - 1) / k
        if (e * d_cand - 1) % k != 0:
            continue
        phi_cand = (e * d_cand - 1) // k

        # Verificar se phi_cand é plausível para n = p*q:
        # phi(n) = n - (p+q) + 1  =>  p+q = n - phi_cand + 1
        # p e q são raízes de x^2 - (p+q)x + n = 0
        soma  = n - phi_cand + 1
        delta = soma * soma - 4 * n
        if delta < 0:
            continue
        sqrt_delta = _isqrt(delta)
        if sqrt_delta * sqrt_delta != delta:
            continue
        p = (soma + sqrt_delta) // 2
        q = (soma - sqrt_delta) // 2
        if p * q == n and miller_rabin(p) and miller_rabin(q):
            if verbose:
                print(f"\n    *** SUCESSO ***")
                print(f"    d encontrado : {d_cand}")
                print(f"    p            : {p}")
                print(f"    q            : {q}")
                print(f"    Verificação  : p * q = {p*q} = n  ✓")
            return d_cand

    if verbose:
        print("    d não encontrado (d provavelmente não é pequeno).")
    return None


def _gerar_rsa_com_d_pequeno(bits_n: int = 256) -> tuple[int, int, int, int, int]:
    """
    Gera parâmetros RSA com d deliberadamente pequeno (< n^{1/4}/3)
    para servir como alvo demonstrativo do ataque de Wiener.

    Retorna (p, q, n, e, d).
    """
    # Estratégia: escolher d pequeno, calcular e = d^{-1} mod phi(n)
    bits_primo = bits_n // 2
    limiar_d   = bits_n // 4 - 2   # d com ~(bits_n/4 - 2) bits: abaixo do limiar

    while True:
        p = gerar_primo(bits_primo)
        q = gerar_primo(bits_primo)
        if p == q:
            continue
        n     = p * q
        phi_n = (p - 1) * (q - 1)

        # Sortear d pequeno (< n^{1/4}/3 ≈ 2^{bits_n/4} / 3)
        limite_d = max(3, (1 << limiar_d) // 3)
        d = secrets.randbelow(limite_d - 2) + 2  # d em [2, limite_d)

        if mdc(d, phi_n) != 1:
            continue
        try:
            e = inverso_modular(d, phi_n)
        except ValueError:
            continue

        if e > 1 and mdc(e, phi_n) == 1:
            return p, q, n, e, d


def demo_wiener() -> None:
    """Demonstra o ataque de Wiener sobre parâmetros RSA vulneráveis."""
    _sep("ATAQUE 2: Wiener – Expoente Privado Pequeno")
    print("""
  Princípio: quando d < n^{1/4}/3, o racional k/d aparece como um
  convergente da fração contínua de e/n. Analisar os convergentes
  permite recuperar d sem fatorar n (Wiener, 1990).

  Defesa: usar d de tamanho comparável a n (gerado via CRT ou
  diretamente). O expoente e=65537 padrão garante que d seja grande.
""")
    print("  Gerando parâmetros RSA com d pequeno (256 bits)... ", end="", flush=True)
    p, q, n, e, d_real = _gerar_rsa_com_d_pequeno(bits_n=256)
    print("concluído.\n")

    print(f"  [Segredo do experimento] d real = {d_real}")
    print(f"  [Segredo do experimento] d (bits) = {d_real.bit_length()}")
    print(f"  n (bits) = {n.bit_length()}")
    print(f"  n^(1/4)/3 ≈ {int(n**0.25)//3}  (d deve ser menor que este valor)\n")
    print("  Executando ataque de Wiener...")

    inicio = time.perf_counter()
    d_enc = wiener_attack(e, n, verbose=True)
    elapsed = time.perf_counter() - inicio

    if d_enc is not None:
        assert d_enc == d_real, "d encontrado difere do real!"
        print(f"\n    Tempo total do ataque : {elapsed*1000:.3f} ms")
        print(f"    d real == d encontrado: {d_enc == d_real}  ✓\n")
    else:
        print("  Ataque não convergiu (d pode não estar no limiar vulnerável).\n")


# ── Ataque 3: Módulo Comum ────────────────────────────────────────────────────

def ataque_modulo_comum(C1: int, C2: int,
                        e1: int, e2: int,
                        n: int,
                        verbose: bool = True) -> Optional[int]:
    """
    Ataque de módulo comum.

    Dado C1 = M^{e1} mod n  e  C2 = M^{e2} mod n,
    com mdc(e1, e2) = 1, recupera M sem conhecer d.

    Usa a identidade de Bezout: s*e1 + t*e2 = 1
    =>  M = C1^s * C2^t mod n

    Retorna M se bem-sucedido, None caso mdc(e1,e2) != 1.
    """
    g, s, t = euclides_estendido(e1, e2)
    if g != 1:
        if verbose:
            print(f"    mdc(e1, e2) = {g} != 1. Ataque inaplicável.")
        return None

    # Tratar expoentes negativos: C^{-k} = (C^{-1})^k mod n
    if s < 0:
        inv_C1 = inverso_modular(C1, n)
        parte1 = exponenciacao_modular(inv_C1, -s, n)
    else:
        parte1 = exponenciacao_modular(C1, s, n)

    if t < 0:
        inv_C2 = inverso_modular(C2, n)
        parte2 = exponenciacao_modular(inv_C2, -t, n)
    else:
        parte2 = exponenciacao_modular(C2, t, n)

    M_rec = parte1 * parte2 % n

    if verbose:
        print(f"    e1 = {e1}, e2 = {e2}")
        print(f"    mdc(e1,e2) = {g}  ✓ (condição para o ataque)")
        print(f"    Coeficientes de Bezout: s={s}, t={t}")
        print(f"    Mensagem recuperada: M = {M_rec}")

    return M_rec


def demo_modulo_comum() -> None:
    """Demonstra o ataque de módulo comum."""
    _sep("ATAQUE 3: Módulo Comum")
    print("""
  Cenário: dois usuários (Alice e Bob) compartilham o mesmo módulo n,
  mas usam expoentes públicos diferentes (e1 e e2 coprimos).
  A mesma mensagem M é enviada para ambos.

  O atacante intercepta C1 e C2 e resolve:
      s*e1 + t*e2 = 1  (identidade de Bezout)
      M = C1^s * C2^t mod n

  Sem precisar fatorar n nem conhecer d1 ou d2.
  Defesa: cada par de usuários DEVE usar módulos distintos.
""")
    # Gerar módulo compartilhado e dois pares de expoentes
    n = 61 * 53   # exemplo do relatório (n=3233), facilmente verificável
    e1, e2 = 17, 23
    phi_n = 60 * 52

    d1 = inverso_modular(e1, phi_n)
    d2 = inverso_modular(e2, phi_n)

    M  = 42    # mensagem secreta
    C1 = exponenciacao_modular(M, e1, n)
    C2 = exponenciacao_modular(M, e2, n)

    print(f"  Parâmetros do experimento:")
    print(f"    n  = {n} (módulo compartilhado)")
    print(f"    M  = {M} (mensagem secreta)")
    print(f"    e1 = {e1}, C1 = M^e1 mod n = {C1}")
    print(f"    e2 = {e2}, C2 = M^e2 mod n = {C2}\n")
    print("  Executando ataque de módulo comum...")

    inicio = time.perf_counter()
    M_rec = ataque_modulo_comum(C1, C2, e1, e2, n, verbose=True)
    elapsed = time.perf_counter() - inicio

    if M_rec == M:
        print(f"\n    *** SUCESSO: M recuperado = {M_rec}  ✓")
        print(f"    Tempo: {elapsed*1000:.4f} ms")
        print(f"    Nem d1 nem d2 foram usados. n não foi fatorado.\n")
    else:
        print(f"    Falha: M recuperado = {M_rec} != {M}\n")


# ── Ataque 4: Hastad (e pequeno, múltiplos destinatários) ────────────────────

def _crt_dois(r1: int, n1: int, r2: int, n2: int) -> tuple[int, int]:
    """CRT para dois módulos coprimos: x ≡ r1 (mod n1), x ≡ r2 (mod n2)."""
    g, s, t = euclides_estendido(n1, n2)
    if g != 1:
        raise ValueError(f"n1={n1} e n2={n2} não são coprimos.")
    N = n1 * n2
    x = (r1 * t * n2 + r2 * s * n1) % N
    return x, N


def crt_geral(residuos: list[int], modulos: list[int]) -> tuple[int, int]:
    """
    Teorema Chinês do Resto para lista de congruências.
    Retorna (x, N) onde N = produto dos módulos.
    """
    x, N = residuos[0], modulos[0]
    for r, m in zip(residuos[1:], modulos[1:]):
        x, N = _crt_dois(x, N, r, m)
    return x % N, N


def hastad_attack(cifrotextos: list[int],
                  modulos: list[int],
                  e: int = 3,
                  verbose: bool = True) -> Optional[int]:
    """
    Ataque de Hastad: recupera M quando M^e é enviado a e destinatários
    sem padding, cada um com módulo diferente.

    Dado C_i = M^e mod n_i (i = 1..e):
      1. CRT => M^e mod (n1*n2*...*ne)
      2. Raiz e-ésima inteira

    Retorna M se a raiz e-ésima for exata, None caso contrário.
    """
    if len(cifrotextos) < e or len(modulos) < e:
        raise ValueError(f"Hastad requer exatamente {e} pares (C_i, n_i).")

    # Passo 1: CRT
    Me_mod_N, N = crt_geral(cifrotextos[:e], modulos[:e])

    if verbose:
        print(f"    e = {e}")
        print(f"    CRT => M^{e} mod N calculado ({N.bit_length()} bits)")

    # Passo 2: raiz e-ésima inteira
    if e == 3:
        M_rec = _raiz_cubica_inteira(Me_mod_N)
    else:
        # Método genérico por Newton-Raphson
        M_rec = int(round(Me_mod_N ** (1 / e)))
        # Ajuste fino
        for delta in range(-3, 4):
            c = M_rec + delta
            if c >= 0 and c ** e == Me_mod_N:
                M_rec = c
                break

    if M_rec ** e == Me_mod_N:
        if verbose:
            print(f"    Raiz {e}a inteira exata encontrada: M = {M_rec}  ✓")
        return M_rec
    else:
        if verbose:
            print(f"    Raiz não é inteira. M pode ser maior que n^(1/e) ou usar padding.")
        return None


def demo_hastad() -> None:
    """Demonstra o ataque de Hastad com e=3 e três destinatários."""
    _sep("ATAQUE 4: Hastad – Expoente Pequeno (e=3)")
    print("""
  Cenário: M é cifrado com e=3 e enviado a 3 destinatários distintos
  (com módulos n1, n2, n3 diferentes), SEM padding.
  O atacante coleta C1, C2, C3 e usa o Teorema Chinês do Resto para
  reconstruir M^3 em Z e depois extrai a raiz cúbica inteira.
  Defesa: SEMPRE usar padding (OAEP) antes de cifrar.
""")
    e = 3
    # Gerar 3 pares de primos distintos (módulos pequenos para velocidade)
    bits_p = 16
    pares = []
    ns    = set()
    while len(pares) < e:
        p = gerar_primo(bits_p)
        q = gerar_primo(bits_p)
        n = p * q
        if n not in ns and p != q:
            pares.append((p, q, n))
            ns.add(n)

    modulos    = [n for (_, _, n) in pares]
    phi_ns     = [(p-1)*(q-1) for (p, q, _) in pares]

    # Mensagem M deve satisfazer M < min(n_i) (e M^3 < n1*n2*n3)
    M_max = min(modulos) - 1
    M = secrets.randbelow(M_max - 10) + 5
    print(f"  Módulos (32 bits cada approx.):")
    for i, (p, q, n) in enumerate(pares, 1):
        print(f"    n{i} = {n}")
    print(f"  Mensagem secreta: M = {M}")

    cifrotextos = [exponenciacao_modular(M, e, n) for n in modulos]
    print(f"  Cifrotextos: {cifrotextos}\n")
    print("  Executando ataque de Hastad...")

    inicio = time.perf_counter()
    M_rec = hastad_attack(cifrotextos, modulos, e=e, verbose=True)
    elapsed = time.perf_counter() - inicio

    if M_rec == M:
        print(f"\n    *** SUCESSO: M recuperado = {M_rec}  ✓")
        print(f"    Tempo: {elapsed*1000:.4f} ms\n")
    else:
        print(f"    Resultado: {M_rec} (esperado: {M})\n")


# ── Resumo de segurança ───────────────────────────────────────────────────────

def resumo_seguranca() -> None:
    _sep("RESUMO: Ataques x Defesas")
    tabela = [
        ("Força Bruta",         "n pequeno (< 60 bits)",      "Módulo n >= 2048 bits"),
        ("Wiener",              "d < n^{1/4}/3",              "Gerar d grande; usar e=65537"),
        ("Módulo Comum",        "n compartilhado, 2 expoentes","Módulo único por entidade"),
        ("Hastad",              "e pequeno, sem padding",      "Padding OAEP obrigatório"),
        ("Bleichenbacher",      "PKCS#1 v1.5, oracle adaptivo","Usar RSA-OAEP (PKCS#1 v2.2)"),
        ("Timing Attack",       "Medição de tempo de d",       "Blinding; implementação em CRT"),
        ("Shor (quântico)",     "Computador quântico tolerante","Migrar para ML-KEM (FIPS 203)"),
    ]
    print()
    linha_fmt = "  {:<25} {:<35} {}"
    print(linha_fmt.format("Ataque", "Pré-condição", "Defesa"))
    print("  " + "─" * 90)
    for ataque, cond, defesa in tabela:
        print(linha_fmt.format(ataque, cond, defesa))
    print()


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print()
    print("=" * 66)
    print("  rsa_attacks.py — Ataques Clássicos ao RSA")
    print("  Atividade Final – Introdução à Cibersegurança")
    print("  SENAI CIMATEC – Pós-graduação em Comunicação Quântica")
    print("=" * 66)
    print()

    demo_forca_bruta()
    _sep()
    demo_wiener()
    _sep()
    demo_modulo_comum()
    _sep()
    demo_hastad()
    _sep()
    resumo_seguranca()

    print("  Todos os ataques concluídos com sucesso.")
    print()
