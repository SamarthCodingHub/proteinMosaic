import streamlit as st

def home_page():
    # Eye-catching Title with Emoji
    st.markdown("<h1 style='text-align: center; color: #FF6F61;'>🧬 Protein Molecule Mosaic 🧬</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center; color: #A7FFEB;'>Explore and analyze protein structures interactively!</h3>", unsafe_allow_html=True)

    # Columns for Objectives & Features
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

    # About Me Section with Social Links
    st.markdown("---")
    st.header("👨‍💻 About Me")
    st.image("https://avatars.githubusercontent.com/u/203984900?v=4", width=120)
    st.markdown("""
    **Samarth Satalinga Kittad**  
    Passionate developer & computer-aided drug discovery enthusiast.  
    I created this app to make protein analysis accessible, interactive, and visually engaging for everyone!

    [![GitHub](https://img.shields.io/badge/GitHub-Profile-informational?style=flat&logo=github)](https://github.com/samarthskittad)

    📧 [samarthkittad8088@gmail.com](mailto:samarthkittad8088@gmail.com)
    """)

    
    # Optional: Add a nice footer
    st.markdown(
        "<hr style='border: 1px solid #A7FFEB;'>"
        "<p style='text-align: center; color: #A7FFEB;'>Made with ❤️ using Streamlit</p>",
        unsafe_allow_html=True
    )

home_page()
