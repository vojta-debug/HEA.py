import streamlit as st
import numpy as np
import pandas as pd

# --- NASTAVENÍ STRÁNKY ---
st.set_page_config(page_title="HEA Material Designer", layout="wide")

# --- DATA ---
ELEMENT_DATA = {
    'Mg': {'r': 1.60, 'Tm': 923, 'Tb': 1363, 'M': 24.305},
    'Sc': {'r': 1.62, 'Tm': 1814, 'Tb': 3103, 'M': 44.956},
    'Ti': {'r': 1.47, 'Tm': 1941, 'Tb': 3560, 'M': 47.867},
    'Zn': {'r': 1.34, 'Tm': 693, 'Tb': 1180, 'M': 65.380}
}

MIXING_ENTHALPY = {
    ('Mg', 'Sc'): 0, ('Mg', 'Ti'): 16, ('Mg', 'Zn'): -4,
    ('Sc', 'Ti'): 0, ('Sc', 'Zn'): -13,
    ('Ti', 'Zn'): -5
}

def calculate_properties(comp, temp_k=None):
    elements = list(comp.keys())
    x = np.array([comp[el] for el in elements])
    if np.sum(x) == 0: return 0, 0, 0, 0
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
    # Pokud je zadaná teplota, počítáme Omegu pro ni, jinak pro Tm_avg
    t_use = temp_k if temp_k else tm_avg
    omega = (t_use * ds_mix) / (abs(dh_mix * 1000) + 1e-12)
    
    return ds_mix, dh_mix, delta, omega

def at_to_grams(comp_at, total_mass_g):
    # Přepočet na hmotnostní zlomky
    molar_masses = {el: ELEMENT_DATA[el]['M'] for el in comp_at}
    total_molar_weight = sum(comp_at[el] * molar_masses[el] for el in comp_at)
    if total_molar_weight == 0: return {el: 0.0 for el in comp_at}
    
    grams = {}
    for el in comp_at:
        wt_frac = (comp_at[el] * molar_masses[el]) / total_molar_weight
        grams[el] = round(wt_frac * total_mass_g, 2)
    return grams

# --- UI ---
st.title("🧪 HEA Material Designer: Mg-Sc-Ti-Zn")

# Část 0: Parametry prostředí
st.sidebar.header("Nastavení procesu")
temp_c = st.sidebar.number_input("Teplota slinování (°C)", value=500)
temp_k = temp_c + 273.15
total_mass = st.sidebar.number_input("Celková navážka vzorku (g)", value=50.0, step=1.0)

# Část 1: Základní (kolegova) slitina
st.header("1. Referenční slitina (Základ)")
col_base = st.columns(3)
with col_base[0]: b_mg = st.number_input("Base Mg (at. %)", value=40.0)
with col_base[1]: b_ti = st.number_input("Base Ti (at. %)", value=30.0)
with col_base[2]: b_zn = st.number_input("Base Zn (at. %)", value=30.0)
base_comp = {'Mg': b_mg, 'Sc': 0.0, 'Ti': b_ti, 'Zn': b_zn}
_, _, b_delta, b_omega = calculate_properties(base_comp, temp_k)
st.info(f"Referenční hodnoty při {temp_c}°C: Delta = {b_delta:.2f}%, Omega = {b_omega:.2f}")

# Část 2: Optimalizace
st.divider()
st.header("2. Optimalizace složení (Top 5 variant)")
st.write("Hledáme: Max Mg, Min Sc (min 10%), Ti a Zn (min 10%).")

if st.button("🚀 Spustit iterační analýzu"):
    results = []
    # Grid search
    for sc in range(10, 40): # Sc limitujeme shora pro úsporu času
        for ti in range(10, 100 - sc - 10 + 1):
            for zn in range(10, 100 - sc - ti + 1):
                mg = 100 - sc - ti - zn
                if mg < 10: continue # Chceme i nějaký hořčík
                
                comp = {'Mg': mg/100, 'Sc': sc/100, 'Ti': ti/100, 'Zn': zn/100}
                _, _, d, o = calculate_properties(comp, temp_k)
                
                if d < 6.6 and o > 1.1:
                    # Skóre: Priorita Mg, penalizace za Sc
                    score = mg - (sc * 5) 
                    results.append({
                        'Mg (at%)': mg, 'Sc (at%)': sc, 'Ti (at%)': ti, 'Zn (at%)': zn,
                        'Delta (%)': round(d, 2), 'Omega': round(o, 2), 'score': score
                    })
    
    if results:
        # Seřazení a výběr top 5
        df = pd.DataFrame(results).sort_values(by='score', ascending=False).head(5)
        st.success("Nalezeno 5 nejlepších variant vyhovujících stabilitě:")
        
        for i, row in df.iterrows():
            with st.expander(f"Varianta {i+1}: Mg{int(row['Mg (at%)'])} Ti{int(row['Ti (at%)'])} Zn{int(row['Zn (at%)'])} Sc{int(row['Sc (at%)'])}"):
                c1, c2 = st.columns(2)
                with c1:
                    st.write("**Fyzikální parametry:**")
                    st.write(f"Delta: {row['Delta (%)']} %")
                    st.write(f"Omega (při {temp_c}°C): {row['Omega']}")
                
                with c2:
                    current_at = {'Mg': row['Mg (at%)'], 'Sc': row['Sc (at%)'], 'Ti': row['Ti (at%)'], 'Zn': row['Zn (at%)']}
                    grams = at_to_grams(current_at, total_mass)
                    st.write(f"**Navážka pro {total_mass} g:**")
                    st.code(f"Mg: {grams['Mg']}g | Sc: {grams['Sc']}g\nTi: {grams['Ti']}g | Zn: {grams['Zn']}g")
    else:
        st.error("Nebyla nalezena žádná kombinace splňující kritéria stability.")

st.sidebar.markdown("---")
st.sidebar.info("💡 **Poznámka k Mg:** Při teplotě slinování nad 1090°C (Tb) je nutné počítat s vysokou tenzí par v uzavřeném systému.")
