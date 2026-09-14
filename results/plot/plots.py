from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
import csv


# ============================================================
# CAMINHOS
# ============================================================

BASE_DIR = Path(__file__).resolve().parent

INPUT_DIR = BASE_DIR / "vtk"

# Salva os resultados na pasta results
OUTPUT_CSV = BASE_DIR.parent / "poiseuille_results.csv"

# Pasta para os gráficos
PLOT_DIR = BASE_DIR / "plots"
PLOT_DIR.mkdir(exist_ok=True)


# ============================================================
# ESCALAS
# ============================================================

# Por enquanto, trabalhar numericamente em lattice units.

RHO_SCALE = 1.0
U_SCALE = 1.0

DX_PHYS = 1.0
DY_PHYS = 1.0

DEPTH = 1.0


# D2Q9
CS2 = 1.0 / 3.0


# None = calcula a vazão média das seções centrais
SECTION_X = None


# ============================================================
# LEITURA DO VTI
# ============================================================

def read_vti_ascii(path):
    """
    Lê DataArray ASCII dos VTI sem depender de VTK/PyVista.
    """

    text = path.read_text()

    extent = re.search(
        r'WholeExtent="([^"]+)"',
        text
    )

    if not extent:
        raise ValueError(
            f"WholeExtent não encontrado em {path}"
        )

    e = list(
        map(
            int,
            extent.group(1).split()
        )
    )

    xmin, xmax, ymin, ymax, zmin, zmax = e

    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    arrays = {}

    for m in re.finditer(
        r'<DataArray([^>]*)>(.*?)</DataArray>',
        text,
        re.S
    ):

        attrs = m.group(1)
        body = m.group(2)

        nm = re.search(
            r'Name="([^"]+)"',
            attrs
        )

        if not nm:
            continue

        name = nm.group(1)

        ncomp = re.search(
            r'NumberOfComponents="(\d+)"',
            attrs
        )

        ncomp = int(ncomp.group(1)) if ncomp else 1

        vals = np.fromstring(
            body,
            sep=' ',
            dtype=np.float64
        )

        if ncomp > 1:

            vals = vals.reshape(
                -1,
                ncomp
            )

        else:

            # VTK: x varia mais rapidamente
            vals = vals.reshape(
                ny,
                nx
            )

        arrays[name] = vals

    return arrays, nx, ny


# ============================================================
# ANÁLISE DE UM VTK
# ============================================================

def analyse(path):

    a, nx, ny = read_vti_ascii(path)

    rho = a['rho']
    ux = a['ux']
    uy = a['uy']


    # ========================================================
    # COORDENADAS
    # ========================================================

    x = np.arange(nx) * DX_PHYS
    y = np.arange(ny) * DY_PHYS


    # ========================================================
    # DOMÍNIO FLUIDO
    # ========================================================

    # Paredes em:
    #
    # y = 0
    # y = NY-1
    #
    # portanto usamos os pontos internos.

    y_slice = slice(
        1,
        ny - 1
    )


    # ========================================================
    # VAZÃO VOLUMÉTRICA
    # ========================================================

    # Q(x) = integral ux(x,y) dy
    #
    # O resultado é um vetor:
    #
    # q_sections[0]   -> Q_in
    # q_sections[-1]  -> Q_out

    q_sections = np.trapezoid(
        ux[y_slice, :] * U_SCALE,
        x=y[y_slice],
        axis=0
    )


    # Entrada
    Q_in = q_sections[0] * DEPTH

    # Saída
    Q_out = q_sections[-1] * DEPTH


    # Região central do canal
    i0 = max(
        1,
        nx // 4
    )

    i1 = min(
        nx - 1,
        3 * nx // 4
    )


    # Vazão média no trecho central
    Q_mean = (
        np.mean(
            q_sections[i0:i1]
        )
        * DEPTH
    )


    # Desvio padrão da vazão
    Q_std = (
        np.std(
            q_sections[i0:i1]
        )
        * DEPTH
    )


    # ========================================================
    # VAZÃO MÁSSICA
    # ========================================================

    # mdot(x) =
    # integral rho(x,y) * ux(x,y) dy

    mdot_sections = np.trapezoid(
        rho[y_slice, :]
        * RHO_SCALE
        * ux[y_slice, :]
        * U_SCALE,

        x=y[y_slice],

        axis=0
    )


    # Entrada
    mdot_in = mdot_sections[0] * DEPTH

    # Saída
    mdot_out = mdot_sections[-1] * DEPTH


    # Média no trecho central
    mdot_mean = (
        np.mean(
            mdot_sections[i0:i1]
        )
        * DEPTH
    )


    # ========================================================
    # ERRO DE CONSERVAÇÃO
    # ========================================================

    if abs(Q_in) > 1e-15:

        Q_error = (
            abs(Q_out - Q_in)
            / abs(Q_in)
        )

    else:

        Q_error = np.nan


    if abs(mdot_in) > 1e-15:

        mdot_error = (
            abs(mdot_out - mdot_in)
            / abs(mdot_in)
        )

    else:

        mdot_error = np.nan


    # ========================================================
    # PRESSÃO
    # ========================================================

    # D2Q9:
    #
    # p = cs² rho

    p_lattice = CS2 * rho

    p = (
        p_lattice
        * RHO_SCALE
        * U_SCALE**2
    )


    p_mean = np.mean(p)


    # Pressão média na entrada
    p_in = np.mean(
        p[y_slice, 0]
    )


    # Pressão média na saída
    p_out = np.mean(
        p[y_slice, -1]
    )


    delta_p = p_in - p_out


    # ========================================================
    # ENERGIA CINÉTICA
    # ========================================================

    kinetic_density = (
        0.5
        * rho
        * RHO_SCALE
        * (
            (ux * U_SCALE)**2
            +
            (uy * U_SCALE)**2
        )
    )


    # Integral no domínio 2D

    E_k_prime = np.trapezoid(

        np.trapezoid(
            kinetic_density,
            x=x,
            axis=1
        ),

        x=y,

        axis=0
    )


    E_k_total = (
        E_k_prime
        * DEPTH
    )


    # ========================================================
    # ENERGIA CINÉTICA ESPECÍFICA
    # ========================================================

    mass_prime = np.trapezoid(

        np.trapezoid(
            rho * RHO_SCALE,
            x=x,
            axis=1
        ),

        x=y,

        axis=0
    )


    if mass_prime != 0:

        e_k_specific = (
            E_k_prime
            / mass_prime
        )

    else:

        e_k_specific = np.nan


    # ========================================================
    # PRESSÃO DINÂMICA
    # ========================================================

    if SECTION_X is None:

        ix = (
            i0 + i1 - 1
        ) // 2

    else:

        ix = int(SECTION_X)


    p_dyn = 0.5 * np.mean(

        rho[y_slice, ix]
        * RHO_SCALE
        *
        (
            (ux[y_slice, ix] * U_SCALE)**2
            +
            (uy[y_slice, ix] * U_SCALE)**2
        )
    )


    # ========================================================
    # PERFIS
    # ========================================================

    y_half = ny // 2


    # Vazão ao longo do canal

    Q_profile = (
        q_sections
        * DEPTH
    )


    # Vazão mássica ao longo do canal

    mdot_profile = (
        mdot_sections
        * DEPTH
    )


    # Pressão ao longo de x no centro do canal

    p_centerline = p[
        y_half,
        :
    ]


    # Energia cinética ao longo de x no centro

    e_k_centerline = kinetic_density[
        y_half,
        :
    ]


    # Perfil de velocidade na entrada

    ux_inlet = ux[
        y_slice,
        0
    ]


    # Perfil de velocidade na saída

    ux_outlet = ux[
        y_slice,
        -1
    ]


    # ========================================================
    # RESULTADOS
    # ========================================================

    return {

        'file': path.name,

        'nx': nx,
        'ny': ny,

        'x': x,
        'y': y,

        # -------------------------
        # Vazão
        # -------------------------

        'Q_in': Q_in,
        'Q_out': Q_out,
        'Q_mean': Q_mean,
        'Q_std': Q_std,
        'Q_error': Q_error,

        # -------------------------
        # Vazão mássica
        # -------------------------

        'mdot_in': mdot_in,
        'mdot_out': mdot_out,
        'mdot_mean': mdot_mean,
        'mdot_error': mdot_error,

        # -------------------------
        # Perfis
        # -------------------------

        'Q_profile': Q_profile,
        'mdot_profile': mdot_profile,

        # -------------------------
        # Pressão
        # -------------------------

        'p_mean': p_mean,
        'p_in': p_in,
        'p_out': p_out,
        'delta_p': delta_p,

        # -------------------------
        # Energia
        # -------------------------

        'E_kinetic_prime_J_m': E_k_prime,
        'E_kinetic_total_J': E_k_total,
        'e_kinetic_specific_J_kg': e_k_specific,
        'p_dynamic_section': p_dyn,

        # -------------------------
        # Perfis centrais
        # -------------------------

        'p_centerline': p_centerline,
        'e_k_centerline': e_k_centerline,

        'ux_inlet': ux_inlet,
        'ux_outlet': ux_outlet,
    }


# ============================================================
# PRESSÃO
# ============================================================

def plot_pressure_evolution(rows):

    plt.figure(
        figsize=(10, 6)
    )

    for data in rows:

        x = data['x']

        pressure = (
            data['p_centerline']
        )

        step = re.search(
            r'(\d+)',
            data['file']
        ).group(1)

        plt.plot(
            x,
            pressure,
            label=f'{step}'
        )

    plt.xlabel('x')
    plt.ylabel('Pressão')

    plt.title(
        'Pressão ao longo do canal'
    )

    plt.grid(True)

    plt.legend(
        title='Passo'
    )

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / 'pressure_profiles.png',
        dpi=300
    )

    plt.show()


# ============================================================
# VAZÃO VOLUMÉTRICA
# ============================================================

def plot_flow_evolution(rows):

    plt.figure(
        figsize=(10, 6)
    )

    for data in rows:

        x = data['x']

        Q = data['Q_profile']

        step = re.search(
            r'(\d+)',
            data['file']
        ).group(1)

        plt.plot(
            x,
            Q,
            label=f'{step}'
        )

    plt.xlabel('x')
    plt.ylabel('Vazão volumétrica')

    plt.title(
        'Vazão volumétrica ao longo do canal'
    )

    plt.grid(True)

    plt.legend(
        title='Passo'
    )

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / 'flow_profiles.png',
        dpi=300
    )

    plt.show()


# ============================================================
# VAZÃO MÁSSICA
# ============================================================

def plot_mass_flow_evolution(rows):

    plt.figure(
        figsize=(10, 6)
    )

    for data in rows:

        x = data['x']

        mdot = data[
            'mdot_profile'
        ]

        step = re.search(
            r'(\d+)',
            data['file']
        ).group(1)

        plt.plot(
            x,
            mdot,
            label=f'{step}'
        )

    plt.xlabel('x')
    plt.ylabel('Vazão mássica')

    plt.title(
        'Vazão mássica ao longo do canal'
    )

    plt.grid(True)

    plt.legend(
        title='Passo'
    )

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / 'mass_flow_profiles.png',
        dpi=300
    )

    plt.show()


def get_step(filename):
    """Extrai o número do passo do nome do VTK."""
    return int(
        re.search(
            r'(\d+)',
            filename
        ).group(1)
    )


def plot_pressure_time(rows):

    steps = np.array([
        get_step(data['file'])
        for data in rows
    ])

    nx = rows[0]['nx']

    x_positions = [
        0,
        nx // 2,
        nx - 1
    ]

    plt.figure(figsize=(10, 6))

    for x_pos in x_positions:

        pressure = np.array([
            data['p_centerline'][x_pos]
            for data in rows
        ])

        plt.plot(
            steps,
            pressure,
            marker='o',
            label=f'x = {x_pos}'
        )

    plt.xlabel('Passo da simulação')
    plt.ylabel('Pressão')
    plt.title('Pressão ao longo do tempo')

    plt.grid(True)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / 'pressure_time.png',
        dpi=300
    )

    plt.show()


def plot_flow_time(rows):

    steps = np.array([
        get_step(data['file'])
        for data in rows
    ])

    nx = rows[0]['nx']

    x_positions = [
        0,
        nx // 2,
        nx - 1
    ]

    plt.figure(figsize=(10, 6))

    for x_pos in x_positions:

        Q = np.array([
            data['Q_profile'][x_pos]
            for data in rows
        ])

        plt.plot(
            steps,
            Q,
            marker='o',
            label=f'x = {x_pos}'
        )

    plt.xlabel('Passo da simulação')
    plt.ylabel('Vazão volumétrica')
    plt.title('Vazão volumétrica ao longo do tempo')

    plt.grid(True)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / 'flow_time.png',
        dpi=300
    )

    plt.show()


def plot_mass_flow_time(rows):

    steps = np.array([
        get_step(data['file'])
        for data in rows
    ])

    nx = rows[0]['nx']

    x_positions = [
        0,
        nx // 2,
        nx - 1
    ]

    plt.figure(figsize=(10, 6))

    for x_pos in x_positions:

        mdot = np.array([
            data['mdot_profile'][x_pos]
            for data in rows
        ])

        plt.plot(
            steps,
            mdot,
            marker='o',
            label=f'x = {x_pos}'
        )

    plt.xlabel('Passo da simulação')
    plt.ylabel('Vazão mássica')
    plt.title('Vazão mássica ao longo do tempo')

    plt.grid(True)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / 'mass_flow_time.png',
        dpi=300
    )

    plt.show()

# ============================================================
# CSV
# ============================================================

def save_results(rows):

    # Somente grandezas escalares.
    # Não colocamos os vetores Q_profile, p_centerline etc.
    # dentro do CSV.

    fields = [

        'file',

        'Q_in',
        'Q_out',
        'Q_mean',
        'Q_std',
        'Q_error',

        'mdot_in',
        'mdot_out',
        'mdot_mean',
        'mdot_error',

        'p_mean',
        'p_in',
        'p_out',
        'delta_p',

        'E_kinetic_prime_J_m',
        'E_kinetic_total_J',
        'e_kinetic_specific_J_kg',
        'p_dynamic_section',
    ]


    with OUTPUT_CSV.open(
        'w',
        newline=''
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields
        )

        writer.writeheader()

        for row in rows:

            writer.writerow({
                key: row[key]
                for key in fields
            })


# ============================================================
# IMPRESSÃO DOS RESULTADOS
# ============================================================

def print_results(rows):

    print(
        f"{'arquivo':18s}"
        f"{'Q_in':>14s}"
        f"{'Q_out':>14s}"
        f"{'mdot_in':>14s}"
        f"{'mdot_out':>14s}"
        f"{'dp':>14s}"
    )

    for r in rows:

        print(
            f"{r['file']:18s}"
            f"{r['Q_in']:14.6e}"
            f"{r['Q_out']:14.6e}"
            f"{r['mdot_in']:14.6e}"
            f"{r['mdot_out']:14.6e}"
            f"{r['delta_p']:14.6e}"
        )




# ============================================================
# MAIN
# ============================================================

def main():

    files = sorted(

        INPUT_DIR.glob('*.vti'),

        key=lambda p:
            int(
                re.search(
                    r'(\d+)',
                    p.stem
                ).group(1)
            )
    )


    if not files:

        raise FileNotFoundError(
            f"Nenhum .vti encontrado em "
            f"{INPUT_DIR.resolve()}"
        )


    # Analisa todos os VTK

    rows = [
        analyse(f)
        for f in files
    ]


    # Salva CSV

    save_results(rows)


    print(
        f"{len(rows)} arquivos analisados."
    )

    print(
        f"Resultados salvos em: "
        f"{OUTPUT_CSV.resolve()}\n"
    )


    # Mostra resultados

    print_results(rows)


    # ========================================================
    # GRÁFICOS
    # ========================================================

    if rows:

        # Pressão p(x) para todos os passos
        plot_pressure_evolution(rows)

        # Vazão Q(x) para todos os passos
        plot_flow_evolution(rows)

        # Vazão mássica mdot(x) para todos os passos
        plot_mass_flow_evolution(rows)

        plot_pressure_time(rows)

        plot_flow_time(rows)

        plot_mass_flow_time(rows)


# ============================================================
# EXECUÇÃO
# ============================================================

if __name__ == '__main__':
    main()