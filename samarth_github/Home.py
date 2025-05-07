import streamlit as st

def home_page():
    # Eye-catching Title with Emoji
    st.markdown("<h1 style='text-align: center; color: #FF6F61;'>🧬 Protein Molecule Mosaic 🧬</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center; color: #A7FFEB;'>Explore and analyze protein structures interactively!</h3>", unsafe_allow_html=True)

    # Animated Welcome (using st.balloons or st.snow for fun)
    st.balloons()

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

    # Interactive Demo Section
    st.markdown("---")
    st.markdown("#### 🎁 Try it now! Upload a sample PDB file to preview:")
    uploaded_file = st.file_uploader("Choose a PDB file", type=["pdb"])
    if uploaded_file is not None:
        st.success("✅ File uploaded! Go to the Model page for analysis.")
        st.write(f"**Filename:** `{uploaded_file.name}`")
        st.info("Tip: On the Model page, you can visualize and analyze your structure!")

    # Fun Expander for More Info
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
    """)

    # Call to Action Button
    st.markdown("---")
    if st.button("🚀 Go to Model Page"):
        st.success("Navigate to the Model page using the sidebar!")

    # Optional: Add a nice footer
    st.markdown(
        "<hr style='border: 1px solid #A7FFEB;'>"
        "<p style='text-align: center; color: #A7FFEB;'>Made with ❤️ using Streamlit</p>",
        unsafe_allow_html=True
    )

home_page()



