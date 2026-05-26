import streamlit as st
import numpy as np

# --- NASTAVENÍ STRÁNKY ---
st.set_page_config(page_title="HEA Kalkulačka & Optimalizátor", layout="wide")

# --- DATA --- z pubchem
ELEMENT_DATA = {
    'Mg': {'r': 1.50, 'Tm': 650},
    'Sc': {'r': 1.60, 'Tm': 1541},
    'Ti': {'r': 1.40, 'Tm': 1668},
    'Zn': {'r': 1.35, 'Tm': 420}
}

MOLAR_MASS = {
    'Mg': 24.305,
    'Sc': 44.956,
    'Ti': 47.867,
    'Zn': 65.380
}
# z open quantum matherials database https://oqmd.org/materials/composition/Ti%20Zn
MIXING_ENTHALPY = {
    ('Mg', 'Sc'): -4.15, ('Mg', 'Ti'): -0.19, ('Mg', 'Zn'): -10.81,
    ('Sc', 'Ti'): -1.54, ('Sc', 'Zn'): -36.57,
    ('Ti', 'Zn'): -21.42
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

def predict_intermetallics(comp_at, delta, omega):
    """Analyzuje přítomné prvky a predikuje riziko vzniku intermetalik."""
    # Databáze pravděpodobných fází na základě binárních diagramů
    INTERMETALLIC_DATABASE = {
        ('Sc', 'Zn'): {
            'enthalpy': -13,
            'phases': ["ScZn", "ScZn₂", "ScZn₁₂", "Sc₃Zn₁₇"],
            'desc': "Velmi silná afinitní vazba. Při nižších teplotách slinování je vznik těchto fází téměř jistý."
        },
        ('Ti', 'Zn'): {
            'enthalpy': -5,
            'phases': ["TiZn", "TiZn₂", "TiZn₃", "Ti₂Zn₃"],
            'desc': "Mírná tendence k tvorbě intermetalik, zejména pokud je v systému lokální přebytek zinku."
        },
        ('Mg', 'Zn'): {
            'enthalpy': -4,
            'phases': ["MgZn", "MgZn₂", "Mg₂Zn₃", "Mg₇Zn₃"],
            'desc': "Nízká teplota tání těchto fází může způsobit vznik křehkých eutektických struktur na hranicích zrn."
        },
        ('Mg', 'Ti'): {
            'enthalpy': 16,
            'phases': ["Nemá (Segregace fází)"],
            'desc': "Extrémně vysoká kladná entalpie. Tyto dva prvky se nesnáší a budou mít tendenci se od sebe oddělovat (likvační mez)."
        },
        ('Mg', 'Sc'): {
            'enthalpy': 0,
            'phases': ["Mg-Sc (Spojité roztoky)"],
            'desc': "Ideální chování, náchylnost k tvorbě intermetalik je minimální."
        },
        ('Sc', 'Ti'): {
            'enthalpy': -1.54,
            'phases': ["Sc-Ti (Omezený tuhý roztok)"],
            'desc': "Prvky jsou strukturálně blízké, riziko tvorby křehkých intermetalických fází je velmi nízké."
        }
    }

    present_elements = [el for el, val in comp_at.items() if val >= 5.0]
    predictions = []
    
    # Prohledáme všechny dvojice prvků, které mají v materiálu zastoupení aspoň 5 %
    for i in range(len(present_elements)):
        for j in range(i + 1, len(present_elements)):
            pair = tuple(sorted((present_elements[i], present_elements[j])))
            if pair in INTERMETALLIC_DATABASE:
                data = INTERMETALLIC_DATABASE[pair]
                # Riziko stoupá, pokud celková stabilita HEA klesá (vysoká delta nebo nízká omega)
                if data['enthalpy'] < -3 and (delta > 5.0 or omega < 1.5):
                    predictions.append({
                        'pár': f"{pair[0]} - {pair[1]}",
                        'fáze': ", ".join(data['phases']),
                        'popis': data['desc'],
                        'entalpie': data['enthalpy']
                    })
    return predictions

# ==========================================
# ČÁST 1: MANUÁLNÍ KALKULAČKA
# ==========================================
st.title("HEA Kalkulačka")
st.write("Zadejte atomární procenta (at. %).")

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
    with res_col1: st.metric(label="Parametr δ (Parametr nesouladu velikostí atomů ↓, ≤6,6)", value=f"{delta:.2f} %")
    with res_col2: st.metric(label="Parametr Ω (Termodynamický vliv entropie vůči entalpii ↑, ≥1,1)", value=f"{omega:.2f}")

# ==========================================
# ČÁST 2: AUTOMATICKÁ OPTIMALIZACE
# ==========================================
st.divider()
st.title("Část 2: Hledání optimální slitiny")

st.subheader("Nastavení optimalizace")
st.write("Vyberte prvek k zafixování a nastavte teplotu slinování. Aplikace dopočítá zbytek pro dosažení nejvyšší stability (Max Omega, Min Delta).")

col_lock1, col_lock2, col_lock3 = st.columns(3)
with col_lock1:
    locked_element = st.selectbox("Prvek k uzamčení:", ["Mg", "Sc", "Ti", "Zn"])
with col_lock2:
    locked_value = st.number_input(f"Hodnota pro {locked_element} (at. %):", min_value=0.0, max_value=97.0, value=40.0, step=1.0)
with col_lock3:
    temp_c = st.number_input("Teplota slinování (°C)", value=500.0, step=50.0)

temp_k = temp_c + 273.15

if st.button("Spustit optimalizaci"):
    with st.spinner("Iteruji přes všechny možné kombinace..."):
        valid_results = []
        
        other_elements = [e for e in ["Mg", "Sc", "Ti", "Zn"] if e != locked_element]
        el1, el2, el3 = other_elements
        
        remainder = int(100 - locked_value)
        
        for v1 in range(10, remainder - 19): 
            for v2 in range(10, remainder - v1 - 9): 
                v3 = remainder - v1 - v2
                if v3 < 10: continue
                
                comp_at = {locked_element: locked_value, el1: v1, el2: v2, el3: v3}
                comp_frac = {k: v/100 for k, v in comp_at.items()}
                
                # Výpočet pro slinovací teplotu i pro základní teplotu (Tm_avg)
                ds, dh, cur_delta, cur_omega_sinter = calculate_hea_properties(comp_frac, temp_k)
                _, _, _, cur_omega_base = calculate_hea_properties(comp_frac, None)
                
                if cur_delta < 6.6 and cur_omega_sinter > 1.1:
                    score = cur_omega_sinter / (cur_delta + 1e-5)
                    
                    valid_results.append({
                        'comp': comp_at,
                        'props': (cur_delta, cur_omega_sinter, cur_omega_base),
                        'score': score
                    })

        valid_results.sort(key=lambda x: x['score'], reverse=True)
        st.session_state['top_5_results'] = valid_results[:5]

# --- ZOBRAZENÍ VÝSLEDKŮ A NAVÁŽKY ---
if 'top_5_results' in st.session_state and st.session_state['top_5_results']:
    st.success("Nalezeny vyhovující kombinace s optimalizovanou stabilitou!")
    st.divider()
    
    st.subheader("Výpočet laboratorní navážky")
    total_mass = st.number_input("Zadejte celkovou navážku vzorku (g):", min_value=0.1, value=5.0, step=1.0)
    
    st.markdown("### Top 5 doporučených složení (seřazeno podle stability)")
    for idx, res in enumerate(st.session_state['top_5_results']):
        comp = res['comp']
        d, o_sinter, o_base = res['props']
        wt = atomic_to_weight(comp)
        grams = calculate_grams(wt, total_mass)
        
        with st.expander(f"Varianta {idx + 1}: Mg {comp['Mg']} | Sc {comp['Sc']} | Ti {comp['Ti']} | Zn {comp['Zn']} (at. %)", expanded=(idx==0)):
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Hořčík (Mg)", f"{comp['Mg']} at. %", f"{wt['Mg']:.1f} wt. %", delta_color="off")
            c2.metric("Skandium (Sc)", f"{comp['Sc']} at. %", f"{wt['Sc']:.1f} wt. %", delta_color="off")
            c3.metric("Titan (Ti)", f"{comp['Ti']} at. %", f"{wt['Ti']:.1f} wt. %", delta_color="off")
            c4.metric("Zinek (Zn)", f"{comp['Zn']} at. %", f"{wt['Zn']:.1f} wt. %", delta_color="off")
            
            st.write(f"**Vypočtené parametry:** Delta = **{d:.2f} %** | Omega (základní) = **{o_base:.2f}** | Omega (při {temp_c} °C) = **{o_sinter:.2f}**")
            
            st.markdown("##### Laboratorní navážka:")
            st.info(f"**Mg:** {grams['Mg']:.2f} g &nbsp;|&nbsp; **Sc:** {grams['Sc']:.2f} g &nbsp;|&nbsp; **Ti:** {grams['Ti']:.2f} g &nbsp;|&nbsp; **Zn:** {grams['Zn']:.2f} g")

elif 'top_5_results' in st.session_state and not st.session_state['top_5_results']:
    st.error("Při tomto uzamčení prvku a zadaných kritériích stability neexistuje žádná vyhovující kombinace. Zkuste hodnotu změnit.")

# Zavolání nové funkce (dosadíš do ní aktuální složení, deltu a omegu)
im_fases = predict_intermetallics({'Mg': c_mg, 'Sc': c_sc, 'Ti': c_ti, 'Zn': c_zn}, delta, omega)

st.divider()
st.title("Pravděpodobnost vzniku případných intermetalik")
if im_fases:
    st.warning("⚠️ **Detekováno riziko sekundárních fází:**")
    for f in im_fases:
        st.write(f"**Subsystém {f['pár']}:** Může tvořit fáze `{f['fáze']}` (ΔH_mix = {f['entalpie']} kJ/mol). *{f['popis']}*")
else:
    st.success("✅ Podle strukturních kritérií by měla slitina zůstat čistým jednofázovým tuhým roztokem.")

# ==========================================
# ČÁST 3: POUŽITÉ VZORCE          
#===========================================
st.divider()
st.title("Použité vzorce")
#fotky
st.image("1.png", caption="Parametr nesouladu velikostí atomů (co nejnižší)")
st.image("2.png", caption="Průměrný atomový poloměr ve směsi")
st.image("3.png", caption="Termodynamický vliv entropie vůči entalpii pro tvorbu tuhého roztoku (co nejvyšší)")
st.image("4.png", caption="Směšovací entropie")
st.image("5.png", caption="Směšovací entalpie")
st.write("Rozdíl v použité hodnotě teploty, viz níže - Omega (základní) vs Omega při ($T_{slinování}$)")
st.image("6.png", caption="Teplota")
st.image("7.png", caption="Teplota při slinování")
