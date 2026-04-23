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

# Molární hmotnosti (g/mol) pro převod at. % -> wt. %
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

# --- FUNKCE PRO VÝPOČET VLASTNOSTÍ ---
def calculate_hea_properties(comp):
    elements = list(comp.keys())
    x = np.array([comp[el] for el in elements])
    
    # Ochrana proti dělení nulou, pokud by všechno bylo 0
    if np.sum(x) == 0:
        return 0, 0, 0, 0

    x = x / np.sum(x) # Normalizace na zlomek 0-1
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

# --- FUNKCE PRO PŘEVOD ATOMÁRNÍCH % NA HMOTNOSTNÍ % ---
def atomic_to_weight(comp_at):
    total_weight = sum(comp_at[el] * MOLAR_MASS[el] for el in comp_at)
    if total_weight == 0:
        return {el: 0.0 for el in comp_at}
    return {el: (comp_at[el] * MOLAR_MASS[el] / total_weight) * 100 for el in comp_at}


# ==========================================
# ČÁST 1: MANUÁLNÍ KALKULAČKA
# ==========================================
st.title("🔬 Část 1: HEA Kalkulačka")
st.write("Zadejte atomární procenta (at. %). Aplikace automaticky dopočítá hmotnostní procenta (wt. %) a stabilitu.")

# Vstupní pole (brackets) místo posuvníků
col1, col2, col3, col4 = st.columns(4)
with col1: c_mg = st.number_input("Mg (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)
with col2: c_sc = st.number_input("Sc (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)
with col3: c_ti = st.number_input("Ti (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)
with col4: c_zn = st.number_input("Zn (at. %)", min_value=0.0, max_value=100.0, value=25.0, step=1.0)

total_at = c_mg + c_sc + c_ti + c_zn

if total_at == 0:
    st.error("Součet atomárních procent nesmí být nula!")
else:
    # Přepočet na zlomky pro výpočet vlastností
    comp_fractions = {'Mg': c_mg / total_at, 'Sc': c_sc / total_at, 'Ti': c_ti / total_at, 'Zn': c_zn / total_at}
    
    # Výpočet hmotnostních procent
    wt_pct = atomic_to_weight({'Mg': c_mg, 'Sc': c_sc, 'Ti': c_ti, 'Zn': c_zn})
    
    # Zobrazení hmotnostních procent přímo pod vstupy
    col1.caption(f"Hmotnostní: **{wt_pct['Mg']:.1f} %**")
    col2.caption(f"Hmotnostní: **{wt_pct['Sc']:.1f} %**")
    col3.caption(f"Hmotnostní: **{wt_pct['Ti']:.1f} %**")
    col4.caption(f"Hmotnostní: **{wt_pct['Zn']:.1f} %**")

    if total_at != 100:
        st.info(f"💡 Součet zadaných atomárních % je {total_at} %. Pro výpočet byly hodnoty automaticky znormovány na 100 %.")

    # Výpočet a zobrazení Delta a Omega
    ds, dh, delta, omega = calculate_hea_properties(comp_fractions)
    
    res_col1, res_col2 = st.columns(2)
    with res_col1:
        st.metric(label="Delta (δ)", value=f"{delta:.2f} %")
    with res_col2:
        st.metric(label="Parametr Omega (Ω)", value=f"{omega:.2f}")

    if delta < 6.6 and omega > 1.1:
        st.success("✅ Predikce: Stabilní pevný roztok (SS).")
    else:
        st.warning("⚠️ Predikce: Pravděpodobně vzniknou intermetalika nebo vícefázová struktura.")


# ==========================================
# ČÁST 2: AUTOMATICKÁ OPTIMALIZACE
# ==========================================
st.divider()
st.title("⚙️ Část 2: Hledání optimální slitiny")
st.write("""
Tato sekce projde tisíce kombinací a najde takovou, která:
1. Má **Delta < 6.6** a **Omega > 1.1** (stabilní pevný roztok).
2. Obsahuje **minimum Skandia** (ale minimálně 10 at. %).
3. Obsahuje **maximum Hořčíku**.
""")

if st.button("🚀 Spustit optimalizaci"):
    with st.spinner("Iteruji přes všechny možné kombinace... (může to trvat pár vteřin)"):
        best_comp = None
        best_score = -float('inf')
        best_props = None
        
        # Iterace po 1 % (krok 1). 
        for sc in range(10, 101, 1): # Sc od 10 do 100
            for mg in range(0, 100 - sc + 1, 1):
                for ti in range(10, 100 - sc - mg + 1, 1):
                    for zn in range(10, 100 - sc - mg - ti + 1, 1):
                    
                    comp = {'Mg': mg/100, 'Sc': sc/100, 'Ti': ti/100, 'Zn': zn/100}
                    ds, dh, cur_delta, cur_omega = calculate_hea_properties(comp)
                    
                    # Podmínky pro stabilní pevný roztok
                    if cur_delta < 6.6 and cur_omega > 1.1:
                        # Skórovací systém: Chceme co nejvíc Mg a co nejméně Sc
                        # Skóre = Mg - (velká penalizace za každý % Sc nad 10)
                        score = mg - (sc * 10) 
                        
                        if score > best_score:
                            best_score = score
                            best_comp = {'Mg': mg, 'Sc': sc, 'Ti': ti, 'Zn': zn}
                            best_props = (cur_delta, cur_omega)

        if best_comp is not None:
            st.success("🎉 Nalezeno optimální složení!")
            
            # Zobrazení výsledků optimalizace
            opt_wt = atomic_to_weight(best_comp)
            
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Hořčík (Mg)", f"{best_comp['Mg']} at. %", f"{opt_wt['Mg']:.1f} wt. %", delta_color="off")
            c2.metric("Skandium (Sc)", f"{best_comp['Sc']} at. %", f"{opt_wt['Sc']:.1f} wt. %", delta_color="off")
            c3.metric("Titan (Ti)", f"{best_comp['Ti']} at. %", f"{opt_wt['Ti']:.1f} wt. %", delta_color="off")
            c4.metric("Zinek (Zn)", f"{best_comp['Zn']} at. %", f"{opt_wt['Zn']:.1f} wt. %", delta_color="off")
            
            st.write(f"**Vypočtené parametry:** Delta = {best_props[0]:.2f} %, Omega = {best_props[1]:.2f}")
        else:
            st.error("❌ Při zadaných kritériích (Sc >= 10%, stabilní roztok) neexistuje žádná vyhovující kombinace. Zkuste zmírnit limity.")
