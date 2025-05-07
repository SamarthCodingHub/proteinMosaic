import random
from streamlit_lottie import st_lottie

def load_lottieurl(url):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

def sidebar_controls():
    # Lottie animation for sidebar (science themed)
    lottie_url = "https://assets10.lottiefiles.com/packages/lf20_3vbOcw.json"
    lottie_json = load_lottieurl(lottie_url)

    fun_facts = [
        "Titin is the largest known protein with nearly 27,000 amino acids!",
        "Hemoglobin carries oxygen in your blood.",
        "Insulin was the first protein to be sequenced.",
        "Some proteins can repair themselves after damage!",
        "Collagen is the most abundant protein in mammals.",
        "Green Fluorescent Protein (GFP) glows under UV light!"
    ]

    with st.sidebar:
        # Welcome animation and fun fact
        if lottie_json:
            st_lottie(lottie_json, height=80, key="sidebar_lottie")
        st.title("Protein Molecule Mosaic")
        st.info("💡 " + random.choice(fun_facts))

        st.markdown("---")

        # Quick navigation
        st.markdown("### 🚀 Quick Navigation")
        st.button("Go to Home Page", on_click=lambda: st.experimental_rerun())
        st.markdown("[What is a PDB file?](https://www.rcsb.org/)", unsafe_allow_html=True)

        st.markdown("---")

        # Theme/Color customizer
        st.markdown("### 🎨 Theme Customizer")
        color = st.color_picker("Pick 3D viewer background", "#ffffff")

        st.markdown("---")

        # Protein fun poll
        st.markdown("### 🗳️ Protein Poll")
        poll = st.radio("Which protein fascinates you most?", ["Hemoglobin", "Insulin", "Collagen", "Titin", "GFP"])
        if st.button("Vote!"):
            st.success(f"You voted for {poll}!")

        st.markdown("---")

        # Mini protein gallery
        st.markdown("### 🖼️ Protein Gallery")
        gallery = {
            "Hemoglobin (1A3N)": "1A3N",
            "Insulin (4INS)": "4INS",
            "Collagen (1CAG)": "1CAG",
            "GFP (1EMA)": "1EMA"
        }
        for label, pid in gallery.items():
            if st.button(f"Load {label}"):
                st.session_state['pdb_id'] = pid
                st.experimental_rerun()

        st.markdown("---")

        # Sidebar mutation simulator preview
        st.markdown("### 🧬 Quick Mutation")
        chain = st.text_input("Chain", "A", max_chars=1, key="sidebar_chain")
        resnum = st.number_input("Residue #", min_value=1, max_value=10000, value=10, key="sidebar_resnum")
        aa = st.selectbox("Mutate to", [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL"
        ], key="sidebar_aa")
        if st.button("Simulate Mutation (Sidebar)"):
            st.success(f"Would mutate Chain {chain} Residue {resnum} to {aa}")

        st.markdown("---")

        # Ligand and rendering controls (from your original)
        st.markdown("**Ligand Display Options**")
        show_ligands = st.checkbox("Highlight Ligands", True)
        render_style = st.selectbox(
            "Rendering Style:",
            ["cartoon", "surface", "sphere"],
            index=0,
            help="Choose molecular representation style"
        )

        st.markdown("---")

        # Feedback/contact
        st.markdown("### 💬 Feedback")
        st.markdown(
            "[Send Feedback](mailto:samarthkittad8088@gmail.com?subject=Protein%20Molecule%20Mosaic%20Feedback)",
            unsafe_allow_html=True
        )

        # Return controls for main app logic
        return {
            'render_style': render_style,
            'show_ligands': show_ligands,
            'color': color
        }
