import streamlit as st
import numpy as np

# --- NASTAVENÍ STRÁNKY ---
st.set_page_config(page_title="HEA Kalkulačka", layout="centered")
st.title("🔬 Výpočet parametrů slitiny Mg-Sc-Ti-Zn")
st.write("Upravte zastoupení prvků (v %) a aplikace okamžitě přepočítá hodnoty Delta a Omega.")

# --- DATA (BEZ KŘEMÍKU) ---
ELEMENT_DATA = {
    'Mg': {'r': 1.60, 'Tm': 923},
    'Sc': {'r': 1.62, 'Tm': 1814},
    'Ti': {'r': 1.47, 'Tm': 1941},
    'Zn': {'r': 1.34, 'Tm': 693}
}

# Entalpie míšení (pouze relevantní dvojice)
MIXING_ENTHALPY = {
    ('Mg', 'Sc'): 0, ('Mg', 'Ti'): 16, ('Mg', 'Zn'): -4,
    ('Sc', 'Ti'): 0, ('Sc', 'Zn'): -13,
    ('Ti', 'Zn'): -5
}

# --- FORMULÁŘ PRO VSTUPY ---
st.subheader("Složení slitiny")
col1, col2 = st.columns(2)

with col1:
    c_mg = st.slider("Hořčík (Mg) %", 0, 100, 25)
    c_sc = st.slider("Skandium (Sc) %", 0, 100, 25)
with col2:
    c_ti = st.slider("Titan (Ti) %", 0, 100, 25)
    c_zn = st.slider("Zinek (Zn) %", 0, 100, 25)

# Celková suma a normalizace
total = c_mg + c_sc + c_ti + c_zn

if total == 0:
    st.error("Součet procent nesmí být nula!")
else:
    # Převod na zlomky (aby součet byl vždy 1.0)
    comp = {
        'Mg': c_mg / total,
        'Sc': c_sc / total,
        'Ti': c_ti / total,
        'Zn': c_zn / total
    }

    # --- VÝPOČET ---
    elements = list(comp.keys())
    x = np.array([comp[el] for el in elements])
    r = np.array([ELEMENT_DATA[el]['r'] for el in elements])
    tm = np.array([ELEMENT_DATA[el]['Tm'] for el in elements])

    # 1. dS_mix (Entropie)
    R = 8.314
    ds_mix = -R * np.sum(x * np.log(x + 1e-12))

    # 2. Delta (Velikostní nesoulad)
    r_avg = np.sum(x * r)
    delta = np.sqrt(np.sum(x * (1 - r / r_avg) ** 2)) * 100

    # 3. dH_mix (Entalpie)
    dh_mix = 0
    for i in range(len(elements)):
        for j in range(i + 1, len(elements)):
            pair = tuple(sorted((elements[i], elements[j])))
            h_ij = MIXING_ENTHALPY.get(pair, 0)
            dh_mix += 4 * h_ij * x[i] * x[j]

    # 4. Omega (Stabilita)
    tm_avg = np.sum(x * tm)
    omega = (tm_avg * ds_mix) / (abs(dh_mix * 1000) + 1e-12)

    # --- ZOBRAZENÍ VÝSLEDKŮ ---
    st.divider()
    if total != 100:
        st.info(f"Poznámka: Součet zadání je {total}%, hodnoty byly automaticky přepočítány na 100%.")

    res_col1, res_col2 = st.columns(2)
    
    with res_col1:
        st.metric(label="Delta (δ)", value=f"{delta:.2f} %")
        st.caption("Limit pro pevné roztoky: < 6.6 %")

    with res_col2:
        st.metric(label="Parametr Omega (Ω)", value=f"{omega:.2f}")
        st.caption("Stabilní roztok při: > 1.1")

    # Interpretace výsledků
    if delta < 6.6 and omega > 1.1:
        st.success("✅ Tato kombinace má vysoký předpoklad pro vznik stabilního pevného roztoku.")
    else:
        st.warning("⚠️ Tato kombinace pravděpodobně vytvoří vícefázovou strukturu nebo intermetalika.")
