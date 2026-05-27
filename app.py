import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

try:
    from Bio import Entrez, SeqIO, Align
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqUtils import gc_fraction
    BIO_OK = True
except ImportError:
    BIO_OK = False

st.set_page_config(
    page_title="KOMOGEN",
    page_icon="🦎",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap');
html, body, [class*="css"] { font-family: 'Inter', sans-serif; }

.hero {
    background: linear-gradient(135deg, #0a0f2e 0%, #2e3b7d 60%, #1a2459 100%);
    border-radius: 16px;
    padding: 2.5rem 2rem 2rem;
    margin-bottom: 1.5rem;
    text-align: center;
    color: white;
}
.hero h1 { font-size: 2.6rem; font-weight: 700; margin: 0.2rem 0; letter-spacing: -1px; }
.hero .sub { font-size: 0.95rem; opacity: 0.8; margin-top: 0.4rem; max-width: 600px; margin-left: auto; margin-right: auto; }
.badge {
    display: inline-block;
    background: rgba(255,255,255,0.12);
    border: 1px solid rgba(255,255,255,0.25);
    border-radius: 99px;
    padding: 3px 12px;
    font-size: 0.75rem;
    margin: 4px 3px 0;
}

.card {
    background: #fff;
    border: 1px solid #d0d8f0;
    border-radius: 12px;
    padding: 1.3rem 1.5rem;
    margin-bottom: 1rem;
    box-shadow: 0 1px 6px rgba(30,50,130,0.06);
}
.card h4 { color: #2e3b7d; margin: 0 0 0.5rem; font-size: 1rem; }

.metric-box {
    background: linear-gradient(135deg, #eef1fb, #f5f7fd);
    border-left: 4px solid #2e3b7d;
    border-radius: 10px;
    padding: 1rem;
    text-align: center;
    height: 100%;
}
.metric-box .val { font-size: 1.9rem; font-weight: 700; color: #2e3b7d; line-height: 1.1; }
.metric-box .lbl { font-size: 0.75rem; color: #555; margin-top: 5px; }

.seq-box {
    font-family: 'Courier New', monospace;
    font-size: 0.78rem;
    background: #f5f7fd;
    border: 1px solid #d0d8f0;
    border-radius: 8px;
    padding: 0.9rem 1rem;
    overflow-x: auto;
    overflow-y: auto;
    white-space: pre-wrap;
    word-break: break-all;
    color: #2e3b7d;
    max-height: 130px;
    line-height: 1.6;
}

.align-box {
    font-family: 'Courier New', monospace;
    font-size: 0.75rem;
    background: #0d1433;
    color: #8fa8d6;
    border-radius: 8px;
    padding: 1rem 1.2rem;
    overflow-x: auto;
    white-space: pre;
    line-height: 1.9;
    max-height: 260px;
    overflow-y: auto;
}

.step-row {
    display: flex;
    align-items: stretch;
    gap: 8px;
    margin: 1rem 0;
}
.step-card {
    flex: 1;
    background: #fff;
    border: 1px solid #d0d8f0;
    border-radius: 10px;
    padding: 1rem;
    text-align: center;
}
.step-card .icon { font-size: 2rem; }
.step-card .title { font-weight: 600; color: #2e3b7d; font-size: 0.9rem; margin: 4px 0; }
.step-card .desc  { font-size: 0.75rem; color: #666; }
.arrow { display: flex; align-items: center; font-size: 1.8rem; color: #a0b0e0; padding-top: 0; }

.section-title { font-size: 1.05rem; font-weight: 600; color: #2e3b7d; margin: 1.2rem 0 0.5rem; }
.hint { font-size: 0.82rem; color: #888; font-style: italic; margin-bottom: 0.6rem; }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero">
  <div style="font-size:3.2rem; margin-bottom:0.2rem;">🦎</div>
  <h1>KomodoSeq</h1>
  <div class="sub">Komputasi Komparatif Sekuens DNA & RNA Gen Peptida Antimikroba<br>pada Komodo sebagai Upaya Konservasi Genetik</div>
  <div style="margin-top:1rem;">
    <span class="badge">🧬 Bioinformatika</span>
    <span class="badge">📘 IF3211 Komputasi Domain Spesifik</span>
    <span class="badge">🌿 Konservasi Genetik</span>
    <span class="badge">🔬 Biopython + Streamlit</span>
  </div>
</div>
""", unsafe_allow_html=True)

if not BIO_OK:
    st.error("❌ Biopython tidak ditemukan. Jalankan: `pip install biopython`")
    st.stop()

with st.sidebar:
    st.markdown("## ⚙️ Konfigurasi")

    with st.expander("📡 NCBI & Accession ID", expanded=True):
        email_ncbi = st.text_input("Email NCBI", value="18223133@std.stei.itb.ac.id",
                                   help="Wajib diisi sesuai kebijakan NCBI Entrez")
        id_komodo  = st.text_input("Accession — Komodo", value="HQ014620.1")
        id_kadal   = st.text_input("Accession — Kadal", value="XM_008113450.3")

    with st.expander("🔬 Parameter Analisis", expanded=True):
        k_mer_size     = st.slider("Ukuran K-mer", 2, 6, 3)
        seq_len_sample = st.slider("Panjang sampel Dogma Sentral (nt)", 30, 200, 90, step=10)
        align_len      = st.slider("Panjang sekuens Alignment (nt)", 50, 300, 100, step=10)

    with st.expander("🔧 Skor Alignment", expanded=False):
        match_s    = st.number_input("Match score",    value= 1.0, step=0.5)
        mismatch_s = st.number_input("Mismatch score", value=-1.0, step=0.5)
        open_gap   = st.number_input("Open gap score", value=-0.5, step=0.1)
        extend_gap = st.number_input("Extend gap",     value=-0.1, step=0.05)

    st.markdown("---")
    run_btn = st.button("🚀 Jalankan Analisis", use_container_width=True, type="primary")
    st.markdown('<p style="font-size:0.75rem;color:#aaa;text-align:center;">Data diunduh langsung dari NCBI</p>', unsafe_allow_html=True)

if "data_loaded" not in st.session_state:
    st.session_state.data_loaded = False

if run_btn:
    with st.spinner("🔄 Mengunduh data dari NCBI..."):
        try:
            Entrez.email = email_ncbi
            h1 = Entrez.efetch(db="nucleotide", id=id_komodo, rettype="fasta", retmode="text")
            st.session_state.komodo_record = SeqIO.read(h1, "fasta"); h1.close()
            h2 = Entrez.efetch(db="nucleotide", id=id_kadal,  rettype="fasta", retmode="text")
            st.session_state.kadal_record  = SeqIO.read(h2, "fasta"); h2.close()
            st.session_state.data_loaded   = True
            st.toast("✅ Data berhasil diunduh!", icon="🎉")
        except Exception as e:
            st.error(f"❌ Gagal mengunduh: {e}")
            st.session_state.data_loaded = False

# Tabs 
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📥 Akuisisi Data",
    "📊 Analisis Genomik",
    "🔄 Dogma Sentral",
    "🔗 Sequence Alignment",
    "ℹ️ Tentang",
])

with tab1:
    st.markdown("### 📥 Akuisisi Data Sekuens dari NCBI")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(f"""
        <div class="card">
        <h4>🦎 Komodo Dragon</h4>
        <p style="margin:0;font-size:0.85rem;"><em>Varanus komodoensis</em></p>
        <p style="margin:4px 0 0;font-size:0.82rem;color:#555;">
            Accession: <code>{id_komodo}</code><br>
            Gen: Beta-Defensin-1 (AMP)<br>
            DB: NCBI Nucleotide
        </p>
        </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown(f"""
        <div class="card">
        <h4>🦎 Kadal Pembanding</h4>
        <p style="margin:0;font-size:0.85rem;"><em>Anolis carolinensis</em></p>
        <p style="margin:4px 0 0;font-size:0.82rem;color:#555;">
            Accession: <code>{id_kadal}</code><br>
            Gen: Beta-Defensin (AMP)<br>
            DB: NCBI Nucleotide
        </p>
        </div>""", unsafe_allow_html=True)

    st.markdown("""
    <div class="card">
    <h4>🧬 Mengapa Gen Beta-Defensin?</h4>
    <p style="font-size:0.88rem;margin:0;color:#333;">
    Beta-defensin adalah peptida antimikroba (AMP) yang menjadi komponen utama sistem imun bawaan reptil.
    Pada Komodo, gen ini berperan dalam ketahanan terhadap bakteri patogen di lingkungan liar.
    Studi komparatif gen ini antara Komodo dan kadal lain memberikan wawasan tentang evolusi imunitas
    dan mendukung strategi konservasi genetik spesies yang terancam punah ini.
    </p>
    </div>""", unsafe_allow_html=True)

    if not st.session_state.data_loaded:
        st.info("👈 Isi parameter di sidebar, lalu klik **Jalankan Analisis**.")
    else:
        kr = st.session_state.komodo_record
        lr = st.session_state.kadal_record

        st.success(f"✅ Data berhasil dimuat — **{kr.id}** dan **{lr.id}**")

        ca, cb, cc = st.columns(3)
        ca.markdown(f'<div class="metric-box"><div class="val">{len(kr.seq):,}</div><div class="lbl">Panjang Sekuens Komodo (bp)</div></div>', unsafe_allow_html=True)
        cb.markdown(f'<div class="metric-box"><div class="val">{len(lr.seq):,}</div><div class="lbl">Panjang Sekuens Kadal (bp)</div></div>', unsafe_allow_html=True)
        cc.markdown(f'<div class="metric-box"><div class="val">{abs(len(kr.seq)-len(lr.seq)):,}</div><div class="lbl">Selisih Panjang (bp)</div></div>', unsafe_allow_html=True)

        st.markdown('<p class="section-title">🧬 Sekuens Komodo (200 nt pertama)</p>', unsafe_allow_html=True)
        st.markdown(f'<div class="seq-box">{str(kr.seq[:200])}</div>', unsafe_allow_html=True)
        st.markdown(f'<p class="hint">ID: {kr.description}</p>', unsafe_allow_html=True)

        st.markdown('<p class="section-title">🧬 Sekuens Kadal Pembanding (200 nt pertama)</p>', unsafe_allow_html=True)
        st.markdown(f'<div class="seq-box">{str(lr.seq[:200])}</div>', unsafe_allow_html=True)
        st.markdown(f'<p class="hint">ID: {lr.description}</p>', unsafe_allow_html=True)

with tab2:
    st.markdown("### 📊 Analisis Karakteristik Genomik")

    if not st.session_state.data_loaded:
        st.info("⚠️ Jalankan analisis terlebih dahulu.")
    else:
        kr = st.session_state.komodo_record
        lr = st.session_state.kadal_record
        gc_k = gc_fraction(kr.seq) * 100
        gc_l = gc_fraction(lr.seq) * 100

        st.markdown('<p class="section-title">🧮 GC Content</p>', unsafe_allow_html=True)
        ca, cb, cc = st.columns(3)
        ca.markdown(f'<div class="metric-box"><div class="val">{gc_k:.1f}%</div><div class="lbl">GC — Komodo</div></div>', unsafe_allow_html=True)
        cb.markdown(f'<div class="metric-box"><div class="val">{gc_l:.1f}%</div><div class="lbl">GC — Kadal</div></div>', unsafe_allow_html=True)
        cc.markdown(f'<div class="metric-box"><div class="val">{abs(gc_k-gc_l):.1f}%</div><div class="lbl">Selisih GC</div></div>', unsafe_allow_html=True)

        fig, ax = plt.subplots(figsize=(6, 3.2))
        x = ["Komodo\n(V. komodoensis)", "Kadal\n(A. carolinensis)"]
        bars = ax.bar(x, [gc_k, gc_l], color=["#2e3b7d","#a0b0e0"],
                      edgecolor="white", linewidth=1.5, width=0.45)
        ax.axhline(50, color="#999", linestyle="--", linewidth=1, alpha=0.6, label="GC = 50%")
        for bar, val in zip(bars, [gc_k, gc_l]):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.4,
                    f"{val:.2f}%", ha="center", va="bottom", fontweight="bold", fontsize=11)
        ax.set_ylabel("GC Content (%)", fontsize=10)
        ax.set_title("Perbandingan GC Content", fontsize=12, fontweight="bold")
        ax.set_ylim(0, max(gc_k, gc_l)+12)
        ax.spines[["top","right"]].set_visible(False)
        ax.set_facecolor("#f5f7fd")
        ax.legend(fontsize=9)
        fig.patch.set_facecolor("white")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close()

        st.markdown(f'<p class="section-title">🔬 Analisis {k_mer_size}-mer (Top 10)</p>', unsafe_allow_html=True)

        def hitung_kmer(seq, k):
            s = str(seq)
            return Counter([s[i:i+k] for i in range(len(s)-k+1)])

        km_k = hitung_kmer(kr.seq, k_mer_size)
        km_l = hitung_kmer(lr.seq, k_mer_size)

        fig2, axes = plt.subplots(1, 2, figsize=(11, 3.5))
        for ax, km, label, color in [
            (axes[0], km_k, "Komodo", "#2e3b7d"),
            (axes[1], km_l, "Kadal", "#a0b0e0"),
        ]:
            top10 = km.most_common(10)
            basa, freq = zip(*top10)
            ax.bar(basa, freq, color=color, edgecolor="white", linewidth=1.2)
            ax.set_title(f"Top 10 {k_mer_size}-mer — {label}", fontsize=10, fontweight="bold")
            ax.set_xlabel(f"{k_mer_size}-mer", fontsize=9)
            ax.set_ylabel("Frekuensi", fontsize=9)
            ax.tick_params(axis="x", rotation=45, labelsize=7)
            ax.spines[["top","right"]].set_visible(False)
            ax.set_facecolor("#f5f7fd")
        fig2.patch.set_facecolor("white")
        fig2.tight_layout()
        st.pyplot(fig2)
        plt.close()

        st.markdown('<p class="section-title">🧩 Komposisi Basa Nitrogen</p>', unsafe_allow_html=True)

        def base_comp(seq):
            s = str(seq).upper()
            n = len(s)
            return {b: s.count(b)/n*100 for b in "ATGC"}

        bc_k = base_comp(kr.seq)
        bc_l = base_comp(lr.seq)
        colors_b = ["#FF716E","#AD82EA","#99C488","#EED682"]

        fig3, axes3 = plt.subplots(1, 2, figsize=(9, 4))
        for ax, bc, title in [(axes3[0], bc_k, "Komodo"), (axes3[1], bc_l, "Kadal")]:
            wedges, texts, autotexts = ax.pie(
                [bc[b] for b in "ATGC"], labels=list("ATGC"),
                autopct="%1.1f%%", colors=colors_b,
                startangle=90, pctdistance=0.75,
                wedgeprops=dict(edgecolor="white", linewidth=2))
            for at in autotexts:
                at.set_fontsize(9); at.set_fontweight("bold")
            ax.set_title(title, fontweight="bold", fontsize=11)
        fig3.suptitle("Komposisi Basa Nitrogen (A / T / G / C)", fontsize=12, fontweight="bold")
        fig3.patch.set_facecolor("white")
        fig3.tight_layout()
        st.pyplot(fig3)
        plt.close()


with tab3:
    st.markdown("### 🔄 Simulasi Dogma Sentral")

    if not st.session_state.data_loaded:
        st.info("⚠️ Jalankan analisis terlebih dahulu.")
    else:
        kr = st.session_state.komodo_record

        st.markdown("""
        <div class="card">
        <h4>📖 Alur Informasi Genetik</h4>
        <p style="font-size:0.87rem;margin:0;color:#333;">
        Dogma sentral biologi molekuler menggambarkan alur informasi genetik:
        <strong>DNA → mRNA (Transkripsi)</strong> → <strong>Protein (Translasi)</strong>.
        Berikut adalah simulasi proses tersebut pada fragmen gen Beta-Defensin Komodo.
        </p>
        </div>""", unsafe_allow_html=True)

        dna     = kr.seq[:seq_len_sample]
        mrna    = dna.transcribe()
        protein = mrna.translate()

        ca, cb, cc = st.columns(3)
        ca.markdown(f'<div class="metric-box"><div class="val">{len(dna)}</div><div class="lbl">DNA Template (nt)</div></div>', unsafe_allow_html=True)
        cb.markdown(f'<div class="metric-box"><div class="val">{len(mrna)}</div><div class="lbl">mRNA (nt)</div></div>', unsafe_allow_html=True)
        cc.markdown(f'<div class="metric-box"><div class="val">{len(protein)}</div><div class="lbl">Rantai Protein (aa)</div></div>', unsafe_allow_html=True)

        st.markdown("""
        <div class="step-row">
          <div class="step-card">
            <div class="icon">🧬</div>
            <div class="title">DNA Template</div>
            <div class="desc">Rantai cetakan asli</div>
          </div>
          <div class="arrow">→</div>
          <div class="step-card">
            <div class="icon">📜</div>
            <div class="title">Transkripsi</div>
            <div class="desc">T diganti U → mRNA</div>
          </div>
          <div class="arrow">→</div>
          <div class="step-card">
            <div class="icon">🔩</div>
            <div class="title">Translasi</div>
            <div class="desc">Kodon → Asam Amino</div>
          </div>
          <div class="arrow">→</div>
          <div class="step-card">
            <div class="icon">⚗️</div>
            <div class="title">Protein</div>
            <div class="desc">Rantai polipeptida</div>
          </div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown('<p class="section-title">1️⃣ DNA Template</p>', unsafe_allow_html=True)
        st.markdown(f'<div class="seq-box">{str(dna)}</div>', unsafe_allow_html=True)

        st.markdown('<p class="section-title">2️⃣ mRNA (hasil Transkripsi)</p>', unsafe_allow_html=True)
        st.markdown(f'<div class="seq-box" style="color:#000000;">{str(mrna)}</div>', unsafe_allow_html=True)

        st.markdown('<p class="section-title">3️⃣ Rantai Protein (hasil Translasi)</p>', unsafe_allow_html=True)
        st.markdown(f'<div class="seq-box" style="color:#000000;">{str(protein)}</div>', unsafe_allow_html=True)

        if "*" in str(protein):
            sp = str(protein).index("*")
            st.success(f"✅ Stop kodon ditemukan di posisi aa ke-**{sp+1}**. Protein fungsional: **{sp} asam amino**.")
        else:
            st.warning("⚠️ Tidak ditemukan stop kodon pada fragmen ini — coba perbesar panjang sampel.")

with tab4:
    st.markdown("### 🔗 Sequence Alignment")

    if not st.session_state.data_loaded:
        st.info("⚠️ Jalankan analisis terlebih dahulu.")
    else:
        kr = st.session_state.komodo_record
        lr = st.session_state.kadal_record

        st.markdown("""
        <div class="card">
        <h4>🔍 Pairwise Global Alignment (Needleman-Wunsch)</h4>
        <p style="font-size:0.87rem;margin:0;color:#333;">
        Alignment global membandingkan sekuens dari ujung ke ujung menggunakan algoritma
        <strong>Needleman-Wunsch</strong> dengan penalti gap affine.
        Hasil alignment menunjukkan tingkat homologi dan konservasi evolutif antar kedua spesies.
        </p>
        </div>""", unsafe_allow_html=True)

        seq_a = kr.seq[:align_len]
        seq_b = lr.seq[:align_len]

        aligner = Align.PairwiseAligner()
        aligner.mode           = "global"
        aligner.match_score    = match_s
        aligner.mismatch_score = mismatch_s
        aligner.open_gap_score = open_gap
        aligner.extend_gap_score = extend_gap

        skor   = aligner.score(seq_a, seq_b)
        persen = (skor / max(len(seq_a), len(seq_b))) * 100

        ca, cb, cc = st.columns(3)
        ca.markdown(f'<div class="metric-box"><div class="val">{skor:.1f}</div><div class="lbl">Skor Alignment</div></div>', unsafe_allow_html=True)
        cb.markdown(f'<div class="metric-box"><div class="val">{persen:.1f}%</div><div class="lbl">Persen Kemiripan</div></div>', unsafe_allow_html=True)
        cc.markdown(f'<div class="metric-box"><div class="val">{align_len}</div><div class="lbl">Panjang Sekuens (nt)</div></div>', unsafe_allow_html=True)

        st.markdown('<p class="section-title">📈 Tingkat Kemiripan</p>', unsafe_allow_html=True)
        bar_color = "#2e3b7d" if persen >= 60 else ("#EED682" if persen >= 40 else "#FF716E")
        fig, ax = plt.subplots(figsize=(8, 1.4))
        ax.barh([""], [100], color="#eef1fb", height=0.5, edgecolor="#d0d8f0", linewidth=1)
        ax.barh([""], [max(persen, 0)], color=bar_color, height=0.5)
        ax.text(min(persen+1.5, 105), 0, f"{persen:.1f}%", va="center",
                fontweight="bold", color=bar_color, fontsize=13)
        ax.set_xlim(0, 112); ax.set_xticks([0,25,50,75,100])
        ax.set_xticklabels(["0%","25%","50%","75%","100%"])
        ax.spines[["top","right","left"]].set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_facecolor("white"); fig.patch.set_facecolor("white")
        fig.tight_layout()
        st.pyplot(fig); plt.close()

        st.markdown('<p class="section-title">📋 Alignment Terbaik (preview)</p>', unsafe_allow_html=True)
        try:
            alignments = aligner.align(seq_a, seq_b)
            best = next(iter(alignments))
            align_str  = str(best)
            lines      = align_str.split("\n")
            preview    = "\n".join(lines[:15]) + ("\n..." if len(lines) > 15 else "")
            st.markdown(f'<div class="align-box">{preview}</div>', unsafe_allow_html=True)
        except Exception as e:
            st.warning(f"Tidak dapat menampilkan preview alignment: {e}")

        st.markdown('<p class="section-title">💡 Interpretasi Biologis</p>', unsafe_allow_html=True)
        if persen >= 70:
            st.success("🟢 **Kemiripan Tinggi** — Homologi signifikan menandakan konservasi evolutif kuat pada gen beta-defensin kedua spesies reptil ini.")
        elif persen >= 40:
            st.warning("🟡 **Kemiripan Sedang** — Kesamaan parsial menunjukkan hubungan evolutif dengan divergensi pada beberapa region sekuens.")
        else:
            st.error("🔴 **Kemiripan Rendah** — Divergensi evolutif besar. Perlu analisis filogenetik lebih lanjut.")

with tab5:
    st.markdown("### ℹ️ Tentang Proyek")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("""
        <div class="card">
        <h4>📋 Informasi Proyek</h4>
        <p style="font-size:0.87rem;color:#333;">
        <b>Mata Kuliah:</b> IF3211 Komputasi Domain Spesifik<br>
        <b>Topik:</b> Sekuens DNA &amp; RNA<br>
        <b>Database:</b> NCBI Nucleotide (via Biopython Entrez)<br><br>
        <b>Judul:</b><br>
        <em>Komputasi Komparatif Sekuens DNA dan RNA Gen Peptida Antimikroba
        pada Komodo sebagai Upaya Konservasi Genetik</em>
        </p>
        </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown("""
        <div class="card">
        <h4>🔧 Stack Teknologi</h4>
        <p style="font-size:0.87rem;color:#333;">
        🐍 <b>Python 3</b> — Bahasa pemrograman utama<br>
        🧬 <b>Biopython</b> — Entrez, SeqIO, Align, SeqUtils<br>
        📊 <b>Matplotlib</b> — Visualisasi data<br>
        🌐 <b>Streamlit</b> — Antarmuka web interaktif<br>
        🗄️ <b>NCBI Entrez API</b> — Sumber data sekuens
        </p>
        </div>""", unsafe_allow_html=True)

    st.markdown("""
    <div class="card">
    <h4>🦎 Latar Belakang Ilmiah</h4>
    <p style="font-size:0.87rem;color:#333;line-height:1.7;">
    Komodo (<em>Varanus komodoensis</em>) adalah kadal terbesar di dunia dan merupakan spesies
    terancam punah (<b>Vulnerable</b>, IUCN). Gen beta-defensin pada Komodo mengkode peptida
    antimikroba (AMP) yang memberikan ketahanan luar biasa terhadap infeksi bakteri di habitatnya.
    Komputasi komparatif sekuens gen ini antara Komodo dan spesies kadal lain membuka wawasan tentang
    evolusi sistem imun bawaan reptil, serta berkontribusi pada upaya konservasi genetik melalui
    pemahaman keunikan genom Komodo yang perlu dilestarikan.
    </p>
    </div>""", unsafe_allow_html=True)

    with st.expander("📦 Cara Menjalankan"):
        st.code("""# Install dependensi
pip install streamlit biopython matplotlib

# Jalankan aplikasi
streamlit run app.py""", language="bash")