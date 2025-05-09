import streamlit as st
from streamlit_lottie import st_lottie
import requests

def load_lottieurl(url):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

def protein_sidebar():
    with st.sidebar:
        st.title("🧬 Protein Molecule Mosaic")

        # New Lottie animation related to protein / molecules
        # Example: Protein folding animation
        lottie_url = "https://assets10.lottiefiles.com/packages/lf20_9yq2yq6f.json"  # Protein folding animation
        lottie_json = load_lottieurl(lottie_url)
        if lottie_json:
            st_lottie(lottie_json, height=150, key="protein_sidebar_lottie")

        # Navigation menu (optional, replace with your own pages)
        page = st.radio("Navigate", ["Home", "Model", "About"], index=1, help="Select app section")

        # Group rendering controls inside an expander
        with st.expander("Rendering Options"):
            style = st.selectbox("Style", ["cartoon", "surface", "sphere"], help="Choose molecular visualization style")
            color_scheme = st.selectbox("Color Scheme", ["spectrum", "chain", "element", "white"], help="Select color scheme")
            highlight_ligands = st.checkbox("Highlight Ligands", value=True, key="sidebar_ligands", help="Toggle ligand highlighting")

        # Status placeholder for future use
        st.markdown("---")
        st.text("Status:")
        st.progress(0)  # Update dynamically in your app logic if needed

        # You can add more controls or info here as needed

        return page, style, color_scheme, highlight_ligands

# Example usage in your main app:
if __name__ == "__main__":
    page, style, color_scheme, highlight_ligands = protein_sidebar()
    st.write(f"Selected page: {page}")
    st.write(f"Style: {style}, Color Scheme: {color_scheme}, Highlight Ligands: {highlight_ligands}")
