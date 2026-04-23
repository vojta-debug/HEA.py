import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# --- NASTAVENÍ STRÁNKY ---
st.set_page_config(page_title="HEA Analýza", layout="centered")
st.title("🔬 Analýza stability lehkých slitin (HEA)")
st.write("Tato aplikace počítá vliv koncentrace vybraného prvku na stabilitu slitiny Mg-Sc-Ti-Zn.")

# --- KONFIGURACE A DATA PRO LEHKÉ SLITINY ---
ELEMENT_DATA = {
    'Mg': {'r': 1.60, 'Tm': 923, 'chi': 1.31},
    'Sc': {'r': 1.62, 'Tm': 1814, 'chi': 1.36},
    'Ti': {'r': 1.47, 'Tm': 1941, 'chi': 1.54},
    'Zn': {'r': 1.34, 'Tm': 693, 'chi': 1.65},
    'Si': {'r': 1.18, 'Tm': 1687, 'chi': 1.90}
}

MIXING_ENTHALPY = {
    ('Mg', 'Sc'): 0, ('Mg', 'Ti'): 16, ('Mg', 'Zn'): -4, ('Mg', 'Si'): -7,
    ('Sc', 'Ti'): 0, ('Sc', 'Zn'): -13, ('Sc', 'Si'): -60,
    ('Ti', 'Zn'): -5, ('Ti', 'Si'): -45,
    ('Zn', 'Si'): 7
}

def calculate_hea_properties(composition):
    elements = list(composition.keys())
    x = np.array([composition[el] for el in elements])
    x = x / np.sum(x)

    r = np.array([ELEMENT_DATA[el]['r'] for el in elements])
    tm = np.array([ELEMENT_DATA[el]['Tm'] for el in elements])

    R = 8.314
    ds_mix = -R * np.sum(x * np.log(x + 1e-12))
    r_avg = np.sum(x * r)
    delta = np.sqrt(np.sum(x * (1 - r / r_avg) ** 2)) * 100

    dh_mix = 0
    for i in range(len(elements)):
        for j in range(i + 1, len(elements)):
            pair = tuple(sorted((elements[i], elements[j])))
            h_ij = MIXING_ENTHALPY.get(pair, 0)
            dh_mix += 4 * h_ij * x[i] * x[j]

    tm_avg = np.sum(x * tm)
    omega = (tm_avg * ds_mix) / (abs(dh_mix * 1000) + 1e-12)

    return ds_mix, dh_mix, delta, omega

# --- GENEROVÁNÍ ANALÝZY A GRAFU ---
def run_sensitivity_analysis(target_element='Si'):
    others = [e for e in ELEMENT_DATA.keys() if e != target_element]
    concentrations = np.linspace(0, 0.25, 50)
    results = {'delta': [], 'omega': [], 'dh_mix': []}

    for c in concentrations:
        remaining_share = (1.0 - c) / len(others)
        comp = {target_element: c}
        for el in others:
            comp[el] = remaining_share

        ds, dh, delta, omega = calculate_hea_properties(comp)
        results['delta'].append(delta)
        results['omega'].append(omega)
        results['dh_mix'].append(dh)

    # Vykreslení grafu
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = 'tab:red'
    ax1.set_xlabel(f'Koncentrace {target_element} (atomární zlomek)')
    ax1.set_ylabel('Delta (velikostní nesoulad %)', color=color)
    ax1.plot(concentrations, results['delta'], color=color, linewidth=2, label='Delta')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.axhline(y=6.6, color='red', linestyle='--', alpha=0.3, label='Limit pro SS (6.6%)')

    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Parametr Omega (stabilita)', color=color)
    ax2.plot(concentrations, results['omega'], color=color, linewidth=2, label='Omega')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.axhline(y=1.1, color='blue', linestyle='--', alpha=0.3, label='Limit pro SS (1.1)')

    plt.title(f'Vliv {target_element} na stabilitu slitiny Mg-Sc-Ti-Zn')
    fig.tight_layout()
    plt.grid(alpha=0.2)
    
    # ZOBRAZENÍ VE STREAMLITU (náhrada za plt.show())
    st.pyplot(fig)

if __name__ == "__main__":
    # Výběr prvku přímo na webu
    target = st.selectbox("Vyber prvek pro analýzu citlivosti:", list(ELEMENT_DATA.keys()))
    run_sensitivity_analysis(target)
