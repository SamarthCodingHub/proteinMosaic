import streamlit as st
from streamlit_lottie import st_lottie
import requests
import random

def load_lottieurl(url):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

def home_page():
    # Load Lottie Animation (science themed)
    lottie_url = "https://assets10.lottiefiles.com/packages/lf20_3vbOcw.json"
    lottie_json = load_lottieurl(lottie_url)
    if lottie_json:
        st_lottie(lottie_json, height=200, key="science_animation")

    # Title and Subtitle with centered styling
    st.markdown("<h1 style='text-align: center; color: #FF6F61;'>🧬 Protein Molecule Mosaic 🧬</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center; color: #A7FFEB;'>Explore and analyze protein structures interactively!</h3>", unsafe_allow_html=True)

    # Rotating Fun Protein Facts
    fun_facts = [
        "Titin is the largest known protein with nearly 27,000 amino acids!",
        "Hemoglobin carries oxygen in your blood.",
        "Insulin was the first protein to be sequenced.",
        "Some proteins can repair themselves after damage!",
        "Green fluorescent protein (GFP) glows under UV light and is widely used as a biological marker.",
        "Spider silk proteins are stronger than steel by weight.",
        "Proteins fold into unique 3D shapes to perform their functions.",
        "Enzymes speed up chemical reactions in living cells."
    ]
    st.info("💡 Did you know? " + random.choice(fun_facts))

    # Quick Start Guide in an expander
    with st.expander("🚀 Quick Start Guide"):
        st.markdown("""
        1. Go to the **Model** page using the sidebar.
        2. Upload your PDB file or enter a PDB ID.
        3. Explore interactive 3D visualization and analysis tools!
        4. Try the **Mutation Simulator** to mutate residues and see the effects instantly!
        5. Or predict protein structures directly from sequences using **ESMFold**.
        """)

    # Objectives & Features in two columns with ESMFold highlighted
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("🎯 Objectives")
        st.markdown("""
        - **Accessible** protein structure visualization  
        - **Upload** or fetch PDB files  
        - **Predict 3D structures from amino acid sequences using ESMFold**  
        - **Classify** ligands & predict active sites  
        - **Ramachandran plot** generation  
        - **Simulate mutations** and analyze structural impact
        """)
    with col2:
        st.subheader("✨ Features")
        st.markdown("""
        - 3D visualization (cartoon, surface, sphere)  
        - Upload or fetch by PDB ID  
        - **Predict protein structures from sequences with ESMFold**  
        - Ligand classification  
        - Active site prediction  
        - Ramachandran plot analysis  
        - **Mutation Simulator**: mutate any residue and instantly view effects  
        - Sidebar controls
        """)

    # Protein Gallery in sidebar with buttons to load examples
    with st.sidebar:
        if lottie_json:
            st_lottie(lottie_json, height=100, key="sidebar_lottie")
        st.markdown("### 🖼️ Protein Gallery")
        if st.button("Load Hemoglobin (1A3N)"):
            st.session_state['pdb_id'] = "1A3N"
        if st.button("Load Insulin (4INS)"):
            st.session_state['pdb_id'] = "4INS"
        st.info("💡 " + random.choice(fun_facts))

    # Educational Expander about PDB files
    st.markdown("---")
    with st.expander("ℹ️ What is a PDB file?"):
        st.write("""
        A PDB file contains 3D structural data of proteins and other biological molecules.  
        You can get sample files from the [RCSB PDB](https://www.rcsb.org/) website.
        """)

    # About Me Section
    st.markdown("---")
    st.header("👨‍💻 About Me")
    st.image("https://avatars.githubusercontent.com/u/203984900?v=4", width=120)
    st.markdown("""
    **Samarth Satalinga Kittad**  
    Passionate developer & computer-aided drug discovery enthusiast.  
    I created this app to make protein analysis accessible, interactive, and visually engaging for everyone!

    [![GitHub](https://img.shields.io/badge/GitHub-Profile-informational?style=flat&logo=github)](https://github.com/samarthskittad)  
    📧 [samarthkittad8088@gmail.com](mailto:samarthkittad8088@gmail.com)  
    [![LinkedIn](https://img.shields.io/badge/LinkedIn-blue?logo=linkedin)](https://www.linkedin.com/in/your-linkedin-profile)
    """)

    # Mentorship & Acknowledgements
    st.markdown("---")
    st.header("🌟 Mentorship & Acknowledgements")
    st.markdown("""
    Special gratitude to **Dr. Kushagra Kashyap**,  
    *Assistant Professor (Bioinformatics), Department of Life Sciences, DES Pune University*.

    Dr. Kashyap’s mentorship has been instrumental in shaping this project. His guidance and expertise brought clarity and depth, fostering learning and innovation.

    [Connect on LinkedIn](https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb/)
    """)

    # Future Improvements
    st.markdown("---")
    st.header("🚀 Future Improvements")
    st.markdown("""
    **Exciting features coming soon:**

    - **Protein-Ligand Docking:**  
      Perform interactive docking and visualize binding poses.

    - **Molecular Dynamics (MD) Simulation Analysis:**  
      Analyze protein flexibility and conformational changes over time.

    Stay tuned for updates! Follow the project on [GitHub](https://github.com/samarthskittad) for the latest developments. 🚀
    """)

    # Feedback Link
    st.markdown("💬 [Send Feedback](mailto:samarthkittad8088@gmail.com?subject=Protein%20Molecule%20Mosaic%20Feedback)")

    # Footer
    st.markdown(
        "<hr style='border: 1px solid #A7FFEB;'>"
        "<p style='text-align: center; color: #A7FFEB;'>Made with ❤️ by Samarth Satalinga Kittad, a passionate bioinformatics student | 2025</p>",
        unsafe_allow_html=True
    )


if __name__ == "__main__":
    home_page()
