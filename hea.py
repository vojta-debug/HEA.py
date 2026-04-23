import streamlit as st
import numpy as np

# --- NASTAVENÍ STRÁNKY ---
st.set_page_config(page_title="HEA Kalkulačka & Optimalizátor", layout="centered")

# --- DATA ---
ELEMENT_DATA = {
    'Mg': {'r': 1.60, 'Tm': 923},
    'Sc': {'r': 1.62, 'Tm': 1814},
    'Ti': {'r': 1.47, 'Tm': 1941},
    'Zn': {'r': 1.34, 'Tm': 693}
}

MOLAR_MASS = {
    'Mg': 24.305,
    'Sc': 44.956,
    'Ti': 47.867,
    'Zn': 65.380
}

MIXING_ENTHALPY = {
    ('Mg', 'Sc'): 0, ('Mg', 'Ti'): 16, ('Mg', 'Zn'): -4,
    ('Sc', 'Ti'): 0, ('Sc', 'Zn'): -13,
    ('Ti', 'Zn'): -5
}

# --- FUNKCE ---
def calculate_hea_properties(comp, temp_k=None):
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
    # Pokud je zadaná teplota slinování, použijeme ji pro Omegu
    t_use = temp_k if temp_k else tm_avg
    omega = (t_use * ds_mix) / (abs(dh_mix * 1000) + 1e-12)

    return ds_mix, dh_mix, delta, omega

def atomic_to_weight(comp_at):
    total_weight = sum(comp_at[el] * MOLAR_MASS[el] for el in comp_at)
    if total_weight == 0: return {el: 0.0 for el in comp_at}
    return {el: (comp_at[el] * MOLAR_MASS[el] / total_weight) * 100 for el in comp_at}

def calculate_grams(comp_wt, total_mass_g):
    return {el: (comp_wt[el] / 100) * total_mass_g for el in comp_wt}


# ==========================================
# ČÁST 1: MANUÁLNÍ KALKULAČKA (Nezměněna)
# ==========================================
st.title("🔬 Část 1: HEA Kalkulačka")
st.write("Zadejte atomární procenta (at. %). Aplikace automaticky dopočítá hmotnostní procenta (wt. %) a stabilitu.")

col1, col2, col3, col4 = st.columns(4)
with col1: c_mg = st.number_input("Mg (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)
with col2: c_sc = st.number_input("Sc (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)
with col3: c_ti = st.number_input("Ti (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)
with col4: c_zn = st.number_input("Zn (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)

total_at = c_mg + c_sc + c_ti + c_zn

if total_at == 0:
    st.error("Součet atomárních procent nesmí být nula!")
else:
    comp_fractions = {'Mg': c_mg / total_at, 'Sc': c_sc / total_at, 'Ti': c_ti / total_at, 'Zn': c_zn / total_at}
    wt_pct = atomic_to_weight({'Mg': c_mg, 'Sc': c_sc, 'Ti': c_ti, 'Zn': c_zn})
    
    col1.caption(f"Hmotnostní: **{wt_pct['Mg']:.1f} %**")
    col2.caption(f"Hmotnostní: **{wt_pct['Sc']:.1f} %**")
    col3.caption(f"Hmotnostní: **{wt_pct['Ti']:.1f} %**")
    col4.caption(f"Hmotnostní: **{wt_pct['Zn']:.1f} %**")

    if total_at != 100:
        st.info(f"💡 Součet zadaných atomárních % je {total_at} %. Pro výpočet byly hodnoty automaticky znormovány na 100 %.")

    ds, dh, delta, omega = calculate_hea_properties(comp_fractions)
    
    res_col1, res_col2 = st.columns(2)
    with res_col1: st.metric(label="Delta (δ)", value=f"{delta:.2f} %")
    with res_col2: st.metric(label="Parametr Omega (Ω)", value=f"{omega:.2f}")


# ==========================================
# ČÁST 2: AUTOMATICKÁ OPTIMALIZACE
# ==========================================
st.divider()
st.title("⚙️ Část 2: Hledání optimální slitiny")

st.subheader("Parametry procesu a referenční slitina")
st.write("Vložte parametry již fungující slitiny (bez Sc) a nastavte teplotu pro výpočet parametru Omega.")

c_ref1, c_ref2, c_ref3, c_ref4 = st.columns(4)
with c_ref1: ref_mg = st.number_input("Ref. Mg (at. %)", value=40.0)
with c_ref2: ref_ti = st.number_input("Ref. Ti (at. %)", value=30.0)
with c_ref3: ref_zn = st.number_input("Ref. Zn (at. %)", value=30.0)
with c_ref4: temp_c = st.number_input("Teplota slinování (°C)", value=500.0, step=50.0)

temp_k = temp_c + 273.15
ref_comp = {'Mg': ref_mg/100, 'Sc': 0, 'Ti': ref_ti/100, 'Zn': ref_zn/100}
_, _, ref_d, ref_o = calculate_hea_properties(ref_comp, temp_k)
st.info(f"**Referenční hodnoty (při {temp_c} °C):** Delta = {ref_d:.2f} %, Omega = {ref_o:.2f}")

st.write("""
**Optimalizační kritéria:**
1. Delta < 6.6 a Omega > 1.1.
2. Minimálně 10 at. % od Sc, Ti i Zn.
3. Nejvyšší podíl Mg, nejnižší podíl Sc.
""")

if st.button("🚀 Spustit optimalizaci"):
    with st.spinner("Iteruji přes všechny možné kombinace..."):
        valid_results = []
        
        for sc in range(10, 71): 
            for ti in range(10, 100 - sc - 10 + 1): 
                for zn in range(10, 100 - sc - ti + 1): 
                    mg = 100 - sc - ti - zn
                    if mg <= 0: continue
                    
                    comp = {'Mg': mg/100, 'Sc': sc/100, 'Ti': ti/100, 'Zn': zn/100}
                    ds, dh, cur_delta, cur_omega = calculate_hea_properties(comp, temp_k)
                    
                    if cur_delta < 6.6 and cur_omega > 1.1:
                        score = mg - (sc * 10) 
                        valid_results.append({
                            'comp': {'Mg': mg, 'Sc': sc, 'Ti': ti, 'Zn': zn},
                            'props': (cur_delta, cur_omega),
                            'score': score
                        })

        # Seřazení výsledků podle skóre a uložení do session_state (aby nezmizely při zadání navážky)
        valid_results.sort(key=lambda x: x['score'], reverse=True)
        st.session_state['top_5_results'] = valid_results[:5]

# --- ZOBRAZENÍ VÝSLEDKŮ A NAVÁŽKY ---
if 'top_5_results' in st.session_state and st.session_state['top_5_results']:
    st.success("🎉 Nalezeny vyhovující kombinace!")
    st.subheader("Výpočet navážky")
    total_mass = st.number_input("Zadejte celkovou navážku vzorku (g):", min_value=0.1, value=50.0, step=1.0)
    
    st.markdown("### Top 5 doporučených složení")
    for idx, res in enumerate(st.session_state['top_5_results']):
        comp = res['comp']
        d, o = res['props']
        wt = atomic_to_weight(comp)
        grams = calculate_grams(wt, total_mass)
        
        with st.expander(f"Varianta {idx + 1}: Mg {comp['Mg']} | Sc {comp['Sc']} | Ti {comp['Ti']} | Zn {comp['Zn']} (at. %)", expanded=(idx==0)):
            # Tabulka vlastností
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Hořčík (Mg)", f"{comp['Mg']} at. %", f"{wt['Mg']:.1f} wt. %", delta_color="off")
            c2.metric("Skandium (Sc)", f"{comp['Sc']} at. %", f"{wt['Sc']:.1f} wt. %", delta_color="off")
            c3.metric("Titan (Ti)", f"{comp['Ti']} at. %", f"{wt['Ti']:.1f} wt. %", delta_color="off")
            c4.metric("Zinek (Zn)", f"{comp['Zn']} at. %", f"{wt['Zn']:.1f} wt. %", delta_color="off")
            
            st.write(f"**Vypočtené parametry:** Delta = {d:.2f} %, Omega (při {temp_c} °C) = {o:.2f}")
            
            # Navážka v gramech
            st.markdown("##### Laboratorní navážka:")
            st.info(f"**Mg:** {grams['Mg']:.2f} g &nbsp;|&nbsp; **Sc:** {grams['Sc']:.2f} g &nbsp;|&nbsp; **Ti:** {grams['Ti']:.2f} g &nbsp;|&nbsp; **Zn:** {grams['Zn']:.2f} g")

elif 'top_5_results' in st.session_state and not st.session_state['top_5_results']:
    st.error("❌ Při zadaných kritériích neexistuje žádná vyhovující kombinace.")
