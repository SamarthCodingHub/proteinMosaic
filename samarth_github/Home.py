import streamlit as st
from streamlit_lottie import st_lottie
import requests

def load_lottieurl(url):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

def home_page():
    # Lottie Animation (science themed)
    lottie_url = "https://assets10.lottiefiles.com/packages/lf20_3vbOcw.json"
    lottie_json = load_lottieurl(lottie_url)
    st_lottie(lottie_json, height=200, key="science_animation")

    # Title and Subtitle
    st.markdown("<h1 style='text-align: center; color: #FF6F61;'>🧬 Protein Molecule Mosaic 🧬</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center; color: #A7FFEB;'>Explore and analyze protein structures interactively!</h3>", unsafe_allow_html=True)

    # Quick Start Guide
    with st.expander("🚀 Quick Start Guide"):
        st.markdown("""
        1. Go to the **Model** page using the sidebar.
        2. Upload your PDB file or enter a PDB ID.
        3. Explore interactive 3D visualization and analysis tools!
        """)

    # Fun Fact
    st.info("💡 Did you know? The largest known protein is Titin, which has nearly 27,000 amino acids!")

    # Objectives & Features in Columns
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("🎯 Objectives")
        st.markdown("""
        - **Accessible** protein structure visualization  
        - **Upload** or fetch PDB files  
        - **Classify** ligands & predict active sites  
        - **Ramachandran plot** generation
        """)

    with col2:
        st.subheader("✨ Features")
        st.markdown("""
        - 3D visualization (cartoon, surface, sphere)  
        - Upload or fetch by PDB ID  
        - Ligand classification  
        - Active site prediction  
        - Ramachandran plot analysis  
        - Sidebar controls
        """)

    # Educational Expander about PDB files
    st.markdown("---")
    with st.expander("ℹ️ What is a PDB file?"):
        st.write("""
        A PDB file contains 3D structural data of proteins and other biological molecules.  
        You can get sample files from the [RCSB PDB](https://www.rcsb.org/) website.
        """)

    # About Me Section with GitHub, Gmail, LinkedIn
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

    # Feedback Link
    st.markdown("💬 [Send Feedback](mailto:samarthkittad8088@gmail.com?subject=Protein%20Molecule%20Mosaic%20Feedback)")

    # Footer
    st.markdown(
        "<hr style='border: 1px solid #A7FFEB;'>"
        "<p style='text-align: center; color: #A7FFEB;'>Made with ❤️ by Samarth Satalinga Kittad | 2025</p>",
        unsafe_allow_html=True
    )

home_page()
