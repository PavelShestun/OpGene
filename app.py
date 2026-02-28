import streamlit as st
import os
import sys
import time
import pandas as pd
import plotly.graph_objects as go
import json
from io import StringIO, BytesIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Добавляем путь к src, чтобы Python видел модули пакета opgene
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

# Импорт внутренних модулей проекта
from src.opgene.factory import OrganismFactory
from src.opgene.objectives.codon_usage import CodonAdaptationObjective
from src.opgene.objectives.structure import GcContentObjective, RnaFoldingObjective, RepeatAvoidanceObjective
from src.opgene.objectives.motifs import MotifAvoidanceObjective
from src.opgene.objectives.codon_pair import CodonPairObjective
from src.opgene.objectives.biosecurity import BiosecurityObjective
from src.opgene.algorithms.genetic import GeneticOptimizer
from src.opgene.utils.reporting import ReportGenerator
from src.opgene.utils.genbank_export import GenBankExporter

# --- Конфигурация страницы ---
st.set_page_config(
    page_title="OpGene Elite v2.1", 
    page_icon="🧬", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Расширенный CSS для профессионального интерфейса ---
st.markdown("""
    <style>
    /* Основной фон и шрифты */
    .stApp { background-color: #f8f9fa; }
    
    /* Стилизация вкладок (Tabs) */
    .stTabs [data-baseweb="tab-list"] { gap: 10px; background-color: transparent; }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        background-color: #e9ecef;
        border-radius: 8px 8px 0 0;
        padding: 0 25px;
        font-weight: 600;
        color: #495057;
        transition: all 0.3s;
    }
    .stTabs [aria-selected="true"] {
        background-color: #2e7d32 !important;
        color: white !important;
        box-shadow: 0 -2px 10px rgba(46, 125, 50, 0.2);
    }
    
    /* Карточки метрик (Metrics) */
    [data-testid="stMetricValue"] { font-size: 1.8rem !important; color: #1b5e20; }
    .stMetric {
        background-color: white;
        padding: 20px !important;
        border-radius: 12px !important;
        border: 1px solid #e0e0e0 !important;
        box-shadow: 0 4px 6px rgba(0,0,0,0.02) !important;
    }
    
    /* Стилизация кнопок */
    .stButton > button {
        width: 100%;
        border-radius: 8px !important;
        height: 3.5em !important;
        font-weight: 600 !important;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        transition: all 0.3s !important;
    }
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }
    
    /* Контейнеры для кода */
    code { color: #d63384 !important; font-family: 'Courier New', monospace; }
    
    /* Боковая панель */
    .css-1647965 { background-color: #ffffff; border-right: 1px solid #eee; }
    </style>
""", unsafe_allow_html=True)

# --- Функции визуализации (Plotly) ---

def get_gc_plot(dna_seq, window=25):
    """Строит график локального GC-состава в скользящем окне"""
    pos, vals = [], []
    for i in range(len(dna_seq) - window + 1):
        sub = dna_seq[i:i+window]
        gc = sum(1 for b in sub if b in "GC") / window
        vals.append(gc)
        pos.append(i)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=pos, y=vals, 
        name="Local GC", 
        line=dict(color='#2ecc71', width=3),
        fill='tozeroy',
        fillcolor='rgba(46, 204, 113, 0.1)'
    ))
    # Линия целевого диапазона
    fig.add_hrect(y0=0.45, y1=0.55, fillcolor="#27ae60", opacity=0.1, line_width=0, annotation_text="Target Range")
    fig.update_layout(
        height=350, 
        margin=dict(l=10, r=10, t=40, b=10), 
        xaxis_title="Position (bp)", 
        yaxis_title="GC Content", 
        yaxis_range=[0, 1],
        hovermode="x unified"
    )
    return fig

def get_speed_plot(dna_seq, codon_usage_obj, window=10):
    """Строит график скорости трансляции (Ribosome Speed Heatmap)"""
    speeds = []
    positions = []
    weights = codon_usage_obj.weights # Веса Relative Adaptiveness (w)
    
    for i in range(0, len(dna_seq) - (window*3), 3):
        chunk = [dna_seq[j:j+3] for j in range(i, i + window*3, 3)]
        avg_w = sum(weights.get(c, 0.1) for c in chunk) / window
        speeds.append(avg_w)
        positions.append(i // 3)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=positions, y=speeds, 
        name="Ribosome Speed", 
        line=dict(color='#3498db', width=3),
        fill='tozeroy',
        fillcolor='rgba(52, 152, 219, 0.1)'
    ))
    # Пометка зоны рампы (замедления в начале)
    fig.add_vrect(x0=0, x1=15, fillcolor="orange", opacity=0.1, line_width=0, annotation_text="Ramp")
    
    fig.update_layout(
        height=350, 
        margin=dict(l=10, r=10, t=40, b=10), 
        xaxis_title="Codon Position", 
        yaxis_title="Relative Speed (CAI-local)", 
        yaxis_range=[0, 1.1],
        hovermode="x unified"
    )
    return fig

# --- Ядро оптимизации (Engine) ---

def run_optimization_process(aa_seq, org_name, tax_id, email, params, threat_db=None, is_batch=False):
    """
    Централизованная функция для запуска оптимизации. 
    Используется и для одиночного режима, и для пакетного.
    """
    # st.status создает красивый выпадающий список этапов работы
    with st.status(f"🧬 Processing Sequence for {org_name}", expanded=not is_batch) as status:
        st.write("🔍 Connecting to NCBI Entrez via DataLoader...")
        factory = OrganismFactory(email)
        
        st.write("📊 Fetching Genomic Data & Calculating CPS (Codon Pair Score)...")
        profile = factory.create(org_name, tax_id)
        
        # Определяем режим: Максимизация или Гармонизация
        mode = "maximize" if "Maximize" in params['opt_mode'] else "harmonize"
        
        st.write(f"🛠️ Configuring {len(params)} Objectives in '{mode}' mode...")
        objectives = [
            CodonAdaptationObjective(profile.codon_usage, profile.codon_table_id, mode=mode, weight=params['w_cai']),
            CodonPairObjective(profile.cps_table, weight=params['w_cpb']),
            GcContentObjective(target_range=profile.ideal_gc_range, weight=params['w_gc']),
            MotifAvoidanceObjective(
                forbidden=profile.forbidden_motifs, 
                suppress_cpg=params['suppress_cpg'], 
                avoid_internal_rbs=params['avoid_rbs'], 
                weight=params['w_mot']
            ),
            RnaFoldingObjective(weight=params['w_rna']),
            RepeatAvoidanceObjective(weight=params['w_rep']),
            BiosecurityObjective(external_db=threat_db, weight=10.0)
        ]
        
        st.write("🧬 Initializing Memetic Genetic Algorithm...")
        optimizer = GeneticOptimizer(
            objectives, 
            profile, 
            pop_size=params['pop_size'], 
            generations=params['gens']
        )
        
        # Визуальный прогресс для одиночного режима
        if not is_batch:
            progress_bar = st.progress(0, text="Evolution in progress...")
            for i in range(10):
                time.sleep(0.05)
                progress_bar.progress((i+1)*10)
        
        dna, results = optimizer.run(aa_seq)
        status.update(label="✅ Success: Optimization Finished", state="complete", expanded=False)
        
    return dna, results, profile

# --- Боковая панель (Sidebar) ---

with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/2907/2907031.png", width=70)
    st.title("OpGene Elite v2.1")
    st.caption("Professional Bio-Design Suite")
    
    st.divider()
    
    # 1. Загрузка базы угроз
    st.markdown("#### 🛡️ Biosecurity Database")
    threat_file = st.file_uploader("Upload custom threat signatures (JSON)", type=["json"])
    loaded_threats = {}
    if threat_file:
        try:
            loaded_threats = json.load(threat_file)
            st.success(f"Loaded {len(loaded_threats)} signatures")
        except Exception as e:
            st.error(f"Invalid JSON: {e}")

    # 2. Выбор организма
    st.markdown("#### 🏥 Host Organism")
    presets = {
        "Escherichia coli K-12": "83333",
        "Bacillus subtilis": "1423",
        "Homo sapiens": "9606",
        "CHO (Hamster)": "10029",
        "S. cerevisiae": "4932",
        "Custom (Manual TaxID)": "custom"
    }
    sel_name = st.selectbox("Select Target Host:", list(presets.keys()))
    if sel_name == "Custom (Manual TaxID)":
        final_org_name = st.text_input("Manual Name", "Mycoplasma")
        final_tax_id = st.text_input("TaxID", "2104")
    else:
        final_org_name = sel_name
        final_tax_id = presets[sel_name]

    # 3. Настройка весов и стратегий
    st.markdown("#### ⚖️ Weights & Strategy")
    params = {
        'w_cai': st.slider("Codon Usage (CAI)", 0.0, 5.0, 1.5),
        'w_cpb': st.slider("Codon Pair Bias (CPS)", 0.0, 5.0, 2.0),
        'w_gc': st.slider("GC Balance", 0.0, 5.0, 1.0),
        'w_rna': st.slider("5' mRNA Stability", 0.0, 5.0, 1.5),
        'w_rep': st.slider("Repeat Avoidance", 0.0, 10.0, 3.0),
        'w_mot': st.slider("Motifs & RBS Avoid", 0.0, 10.0, 5.0),
        'pop_size': st.select_slider("Population Size", [10, 20, 50, 100], 50),
        'gens': st.select_slider("Number of Generations", [5, 10, 20, 50], 20),
        'opt_mode': st.radio("Strategy Mode:", ["Maximize Efficiency", "Harmonize (Natural)"]),
        'avoid_rbs': st.checkbox("Avoid Internal RBS", True),
        'suppress_cpg': st.checkbox("Suppress CpG (Eukaryotic)", False)
    }
    
    st.divider()
    user_email = st.text_input("NCBI Entrez Email", "your@lab.com")
    st.caption("OpGene Engine v2.1.0 | SOTA Grade")

# --- Основной интерфейс приложения ---

st.title("🧬 OpGene Elite: SOTA Codon Optimizer")
st.caption("Industrial-grade protein expression optimization & gene synthesis readiness check")
st.markdown("---")

# Создание вкладок для одиночного и пакетного режима
t1, t2 = st.tabs(["🎯 Single Protein Optimization", "📦 Batch Processing (FASTA)"])

# --- ВКЛАДКА 1: Одиночная последовательность ---
with t1:
    col_input, col_info = st.columns([2, 1])
    
    with col_input:
        aa_input = st.text_area(
            "Input Amino Acid Sequence:", 
            height=180, 
            placeholder="Paste protein sequence (e.g., MSKGEEL...)",
            value="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRY"
        )
    
    with col_info:
        st.info("""
        **Bio-Logic Info:**
        - **Maximize**: Лучший выбор для производства простых белков.
        - **Harmonize**: Рекомендуется для сложных человеческих белков в бактериях (предотвращает тельца включения).
        - **CPS**: Оптимизирует пары кодонов для предотвращения заторов рибосом.
        """)
        run_btn = st.button("🚀 Start Evolutionary Optimization", type="primary", use_container_width=True)

    if run_btn:
        # Предварительная очистка последовательности
        aa_clean = "".join(aa_input.split()).upper()
        
        if not aa_clean:
            st.error("❌ Error: Amino acid sequence cannot be empty.")
        else:
            try:
                # Запуск движка
                dna, results, profile = run_optimization_process(
                    aa_clean, final_org_name, final_tax_id, user_email, params, loaded_threats
                )
                
                # Проверка на алерты биобезопасности
                if "THREAT" in str(results['metrics']):
                    st.warning("⚠️ **BIOSECURITY ALERT:** Potential threat signatures or antibiotic markers detected in DNA!")

                st.success("✨ Optimization successfully completed!")
                
                # Панель ключевых метрик (KPIs)
                st.markdown("### 📈 Performance Metrics")
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Total Fitness", f"{results['fitness']:.2f}")
                
                # Парсинг строковых метрик для отображения в карточках
                cai_data = results['metrics'].get('CodonAdaptationObjective', '0:0').split(':')
                m2.metric("CAI / MSE", cai_data[-1].strip() if len(cai_data)>1 else "N/A")
                
                gc_data = results['metrics'].get('GcContentObjective', '0:0').split(':')
                m3.metric("GC Content", gc_data[-1].strip() if len(gc_data)>1 else "N/A")
                
                cps_data = results['metrics'].get('CodonPairObjective', '0:0').split(':')
                m4.metric("Avg CPS (Speed)", cps_data[-1].strip() if len(cps_data)>1 else "N/A")

                # Секция визуализации и данных
                res_tab_v, res_tab_d = st.tabs(["📊 Visual Analysis", "📑 Sequence & Export"])
                
                with res_tab_v:
                    # Графики в две колонки
                    g_col1, g_col2 = st.columns(2)
                    with g_col1:
                        st.plotly_chart(get_gc_plot(dna), use_container_width=True)
                    with g_col2:
                        # Получаем веса через временный объект для графика
                        cai_obj_temp = CodonAdaptationObjective(profile.codon_usage, profile.codon_table_id)
                        st.plotly_chart(get_speed_plot(dna, cai_obj_temp), use_container_width=True)
                    
                    st.subheader("Objective Breakdown")
                    st.table(pd.DataFrame.from_dict(results['metrics'], orient='index', columns=['Status']))
                
                with res_tab_d:
                    st.subheader("Optimized DNA Sequence (CDS)")
                    st.code(dna, language="text")
                    st.caption(f"DNA Length: {len(dna)} bp | Calculated based on standard genetic code.")
                    
                    st.divider()
                    st.subheader("📥 Export Results")
                    dl_col1, dl_col2, dl_col3 = st.columns(3)
                    
                    with dl_col1:
                        fasta_data = f">OpGene_{final_org_name}_v2.1\n{dna}"
                        st.download_button("📥 Download FASTA", fasta_data, f"optimized_{final_org_name}.fasta", use_container_width=True)
                    
                    with dl_col2:
                        # Генерация GenBank файла
                        gb_data = GenBankExporter.create_record(dna, final_org_name, aa_clean, profile)
                        st.download_button("📥 Download GenBank (.gb)", gb_data, f"optimized_{final_org_name}.gb", use_container_width=True)
                        
                    with dl_col3:
                        # Генерация PDF отчета
                        pdf_rep = ReportGenerator.create_pdf(final_org_name, aa_clean, dna, results['metrics'], results['fitness'])
                        st.download_button("📄 Download PDF Report", data=bytes(pdf_rep), file_name=f"Report_{final_org_name}.pdf", use_container_width=True)

            except Exception as e:
                st.error(f"Critical Error during optimization: {str(e)}")
                st.exception(e)

# --- ВКЛАДКА 2: Пакетная обработка (Batch Upload) ---
with t2:
    st.markdown("### 📦 Bulk Sequence Optimization")
    st.write("Upload a multi-FASTA file. OpGene will optimize all sequences sequentially using current sidebar settings.")
    
    batch_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "txt"])
    
    if batch_file:
        if st.button("🛠️ Execute Batch Process", type="primary"):
            # Чтение FASTA записей
            stringio = StringIO(batch_file.getvalue().decode("utf-8"))
            records = list(SeqIO.parse(stringio, "fasta"))
            
            if not records:
                st.error("No valid sequences found in the uploaded file.")
            else:
                st.info(f"Detected {len(records)} sequences. Starting processing...")
                batch_final_results = []
                main_progress = st.progress(0)
                
                for i, record in enumerate(records):
                    # Запуск ядра в пакетном режиме (is_batch=True скрывает лишний UI)
                    dna_b, res_b, _ = run_optimization_process(
                        str(record.seq).upper(), final_org_name, final_tax_id, user_email, params, loaded_threats, is_batch=True
                    )
                    
                    batch_final_results.append({
                        "ID": record.id,
                        "AA Length": len(record.seq),
                        "Final Fitness": f"{res_b['fitness']:.2f}",
                        "Metrics summary": str(res_b['metrics']),
                        "Optimized DNA": dna_b
                    })
                    # Обновление прогресс-бара
                    main_progress.progress((i + 1) / len(records))
                
                st.success("✅ All sequences in the batch have been processed.")
                
                # Отображение таблицы результатов (скрываем длинную колонку ДНК для компактности)
                df_results = pd.DataFrame(batch_final_results)
                st.dataframe(df_results.drop(columns=["Optimized DNA"]), use_container_width=True)
                
                # Кнопка скачивания итоговой таблицы
                csv_out = df_results.to_csv(index=False).encode('utf-8')
                st.download_button(
                    "📥 Download All Results (CSV)", 
                    csv_out, 
                    "opgene_batch_results.csv", 
                    "text/csv", 
                    use_container_width=True
                )