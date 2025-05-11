import streamlit as st
import requests
from Bio.PDB import PDBParser, PPBuilder, PDBIO, Superimposer
from Bio.PDB.PDBExceptions import PDBException
from io import StringIO
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import py3Dmol
import biotite.structure.io as bsio
import random

# ----------------------
# Protein Facts
# ----------------------
PROTEIN_FACTS = [
    "The largest known protein is titin, with over 38,000 amino acids.",
    "Hemoglobin, a protein in your blood, carries oxygen from your lungs to your tissues.",
    "Enzymes are proteins that speed up chemical reactions in living cells.",
    "Green fluorescent protein (GFP) glows under UV light and is widely used as a biological marker.",
    "Collagen is the most abundant protein in mammals and provides structure to skin and bones.",
    "Spider silk proteins are stronger than steel by weight.",
    "Insulin, a small protein, regulates blood sugar levels.",
    "Proteins fold into unique 3D shapes to perform their functions.",
    "Proteins can act as molecular machines, like ATP synthase.",
    "The sequence of amino acids determines a protein's structure and function."
]

# ----------------------
# Helper Functions
# ----------------------
@st.cache_data
def fetch_pdb_data(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.text
    except Exception as e:
        st.error(f"Error fetching PDB data: {str(e)}")
        return None

def esmfold_predict_structure(sequence):
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    if response.status_code != 200:
        st.error("ESMFold API error: " + response.text)
        return None, None
    pdb_string = response.content.decode('utf-8')
    with open('predicted_tmp.pdb', 'w') as f:
        f.write(pdb_string)
    struct = bsio.load_structure('predicted_tmp.pdb', extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)
    return pdb_string, b_value

def classify_ligand(residue):
    resname = residue.get_resname().strip()
    if len(resname) <= 2:
        return 'ion'
    elif any(atom.name in ['OXT', 'ND1', 'NE2'] for atom in residue):
        return 'polydentate'
    return 'monodentate'

def extract_ligands(pdb_data):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("temp", StringIO(pdb_data))
    ligands = {'ion': [], 'monodentate': [], 'polydentate': []}
    for residue in structure.get_residues():
        if residue.id[0] != ' ':
            ligand_type = classify_ligand(residue)
            if ligand_type == 'ion':
                ligands['ion'].append(residue.get_resname())
            else:
                ligands[ligand_type].append({
                    'resname': residue.get_resname(),
                    'chain': residue.parent.id,
                    'resnum': residue.id[1],
                    'type': ligand_type
                })
    return ligands

def predict_active_sites(pdb_data):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("temp", StringIO(pdb_data))
    catalytic_residues = ['HIS', 'ASP', 'GLU', 'SER', 'CYS', 'LYS', 'TYR', 'ARG']
    active_sites = []
    ligand_atoms = []
    for residue in structure.get_residues():
        if residue.id[0] != ' ':
            ligand_atoms.extend(list(residue.get_atoms()))
    for residue in structure.get_residues():
        if residue.id[0] == ' ' and residue.get_resname() in catalytic_residues:
            res_atoms = list(residue.get_atoms())
            for res_atom in res_atoms:
                for lig_atom in ligand_atoms:
                    distance = res_atom - lig_atom
                    if distance < 3.0:
                        active_sites.append({
                            'resname': residue.get_resname(),
                            'chain': residue.parent.id,
                            'resnum': residue.id[1],
                            'distance': f"{distance:.2f} Å"
                        })
                        break
                else:
                    continue
                break
    return active_sites

def visualize_ligand_counts(ligands):
    labels = list(ligands.keys())
    counts = [len(ligands[ligand_type]) for ligand_type in labels]
    fig = go.Figure(data=[
        go.Bar(name='Ligand Counts', x=labels, y=counts)
    ])
    fig.update_layout(title='Ligand Type Counts',
                      xaxis_title='Ligand Type',
                      yaxis_title='Count')
    return fig

def get_phi_psi_angles(pdb_string):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", StringIO(pdb_string))
    phi_psi = []
    for model in structure:
        for chain in model:
            ppb = PPBuilder()
            for pp in ppb.build_peptides(chain):
                angles = pp.get_phi_psi_list()
                for phi, psi in angles:
                    if phi is not None and psi is not None:
                        phi_psi.append((phi * 180.0 / 3.14159, psi * 180.0 / 3.14159))
    return phi_psi

def plot_ramachandran(phi_psi):
    fig, ax = plt.subplots(figsize=(5, 5))
    if phi_psi:
        phi, psi = zip(*phi_psi)
        ax.scatter(phi, psi, s=10)
        ax.set_xlabel("Phi")
        ax.set_ylabel("Psi")
        ax.set_title("Ramachandran Plot")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.grid(True)
    return fig

def ramachandran_region_analysis(phi_psi_list):
    favored = 0
    allowed = 0
    outlier = 0
    total = len(phi_psi_list)
    for phi, psi in phi_psi_list:
        if -160 <= phi <= -40 and -80 <= psi <= -20:
            favored += 1
        elif -180 <= phi <= -40 and 90 <= psi <= 180:
            favored += 1
        elif (-180 <= phi <= -20 and -180 <= psi <= 180):
            allowed += 1
        else:
            outlier += 1
    allowed = allowed - favored
    outlier = total - favored - allowed
    return {
        "favored": 100 * favored / total if total else 0,
        "allowed": 100 * allowed / total if total else 0,
        "outlier": 100 * outlier / total if total else 0,
        "total": total
    }

def show_3d_structure(pdb_data, style='cartoon', highlight_ligands=True, width=800, height=500):
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, 'pdb')
    if style == 'cartoon':
        view.setStyle({'cartoon': {'color': 'spectrum'}})
    elif style == 'surface':
        view.setStyle({'cartoon': {'color': 'white'}})
        view.addSurface(py3Dmol.SAS, {'opacity': 0.7})
    elif style == 'sphere':
        view.setStyle({'sphere': {'colorscheme': 'Jmol'}})
    if highlight_ligands:
        view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
    view.zoomTo()
    view.setBackgroundColor('white')
    st.components.v1.html(view._make_html(), height=height, width=width)

def mutate_residue_in_pdb(pdb_data, chain_id, resnum, new_resname):
    lines = pdb_data.splitlines()
    mutated_lines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            line_chain = line[21]
            try:
                line_resnum = int(line[22:26])
            except ValueError:
                mutated_lines.append(line)
                continue
            if line_chain == chain_id and line_resnum == resnum:
                new_line = line[:17] + new_resname.ljust(3) + line[20:]
                mutated_lines.append(new_line)
            else:
                mutated_lines.append(line)
        else:
            mutated_lines.append(line)
    return "\n".join(mutated_lines)

def superimpose_structures(pdb_data1, pdb_data2):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure1 = parser.get_structure("structure1", StringIO(pdb_data1))
    structure2 = parser.get_structure("structure2", StringIO(pdb_data2))

    def get_ca_atoms_with_id(structure):
        atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() != "HOH" and 'CA' in residue:
                        atoms.append((chain.id, residue.id[1], residue['CA']))
        return atoms

    atoms1_with_id = get_ca_atoms_with_id(structure1)
    atoms2_with_id = get_ca_atoms_with_id(structure2)

    common_atoms1 = []
    common_atoms2 = []
    matched_residues = {}

    for chain1, resnum1, atom1 in atoms1_with_id:
        for chain2, resnum2, atom2 in atoms2_with_id:
            if chain1 == chain2 and resnum1 == resnum2:
                common_atoms1.append(atom1)
                common_atoms2.append(atom2)
                matched_residues[(chain1, resnum1)] = True
                break

    if not common_atoms1 or not common_atoms2 or len(common_atoms1) < 3:
        st.error("Not enough common C-alpha atoms found for reliable superimposition.")
        return None, None

    if len(common_atoms1) != len(common_atoms2):
        st.warning(f"The number of common residues found is {len(common_atoms1)}, which is less than the total residues in each structure. The alignment might not be optimal.")

    superimposer = Superimposer()
    superimposer.set_atoms(common_atoms1, common_atoms2)
    try:
        superimposer.apply(structure2.get_atoms())
    except PDBException as e:
        st.error(f"Superimposition failed: {e}")
        return None, None

    aligned_structure_io = StringIO()
    io.set_structure(structure2)
    io.save(aligned_structure_io)
    aligned_structure_str = aligned_structure_io.getvalue()

    return aligned_structure_str, superimposer.rms

# ----------------------
# UI Components
# ----------------------

def sidebar_controls():
    with st.sidebar:
        st.image("https://media.istockphoto.com/id/1390037416/photo/chain-of-amino-acid-or-bio-molecules-called-protein-3d-illustration.jpg?s=612x612&w=0&k=20&c=xSkGolb7TDjqibvINrQYJ_rqrh4RIIzKIj3iMj4bZqI=", width=400)
        st.title("Protein Molecule Mosaic")
        analysis_type = st.radio(
            "Analysis Mode:",
            ["Single Structure", "Structural Comparison"],
            help="Choose the type of analysis to perform."
        )
        render_style = st.selectbox(
            "Rendering Style:",
            ["cartoon", "surface", "sphere"],
            index=0,
            help="Choose molecular representation style"
        )
        st.markdown("---")
        st.markdown("**Ligand Display Options**")
        show_ligands = st.checkbox("Highlight Ligands", True, key="sidebar_ligand_checkbox")
        return {
            'analysis_type': analysis_type,
            'render_style': render_style,
            'show_ligands': show_ligands,
        }

# ----------------------
# Main App Logic
# ----------------------

def main():
    st.set_page_config(
        page_title="Protein Molecule Mosaic",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    # --- INTERACTIVE FACTS & ESMFold HIGHLIGHT ---
    st.markdown("## 🧬 Welcome to Protein Molecule Mosaic!")
    st.info(
        "✨ **Did you know?** " + random.choice(PROTEIN_FACTS),
        icon="💡"
    )
    st.success(
        "🚀 **You can predict protein 3D structure directly from sequence using ESMFold!**\n"
        "Select **'Predict from Sequence (ESMFold)'** below, paste your sequence, and click 'Predict Structure with ESMFold'.",
        icon="🧠"
    )
    st.markdown("---")

    controls = sidebar_controls()
    col1, col2 = st.columns([3, 1])

    with col1:
        st.header("Protein Palette")

        tab1, tab2 = st.tabs(["Single Structure Analysis", "Structural Comparison"])

        with tab1:
            st.subheader("Analyze Single Structure")
            input_mode_single = st.radio(
                "Choose input mode:",
                ["Upload PDB", "Fetch by PDB ID", "Predict from Sequence (ESMFold)"],
                key="single_input_mode"
            )
            pdb_data_single = None
            source_single = None
            plddt = None

            if input_mode_single == "Upload PDB":
                uploaded_pdb = st.file_uploader("Upload a PDB file", type=["pdb"], key="single_upload")
                if uploaded_pdb is not None:
                    pdb_data_single = uploaded_pdb.read().decode("utf-8")
                    source_single = "upload"
                    st.success("PDB file uploaded and loaded.")
            elif input_mode_single == "Fetch by PDB ID":
                pdb_id = st.text_input("Enter PDB ID:", key="single_pdbid").upper()
                if pdb_id:
                    pdb_data_single = fetch_pdb_data(pdb_id)
                    source_single = "pdbid"
                    if pdb_data_single:
                        st.success(f"PDB ID {pdb_id} loaded from RCSB.")
            elif input_mode_single == "Predict from Sequence (ESMFold)":
                example_seq = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
                sequence = st.text_area("Paste your protein sequence (1-letter code):", value=example_seq, height=200)
                if st.button("Predict Structure with ESMFold"):
                    with st.spinner("Predicting structure..."):
                        pdb_data_single, plddt = esmfold_predict_structure(sequence)
                        source_single = "esmfold"
                    if pdb_data_single:
                        st.success("Structure predicted with ESMFold!")
                        st.info(f"Average pLDDT: {plddt}")
                        st.download_button(
                            label="Download Predicted PDB",
                            data=pdb_data_single,
                            file_name='esmfold_predicted.pdb',
                            mime='text/plain',
                        )
                    else:
                        st.error("Prediction failed.")

            if pdb_data_single:
                st.subheader("3D Structure Viewer")
                show_3d_structure(
                    pdb_data_single,
                    style=controls['render_style'],
                    highlight_ligands=controls['show_ligands']
                )
                with st.expander("Ramachandran Plot"):
                    phi_psi = get_phi_psi_angles(pdb_data_single)
                    if phi_psi:
                        fig = plot_ramachandran(phi_psi)
                        st.pyplot(fig)
                        ramachandran_analysis = ramachandran_region_analysis(phi_psi)
                        st.write("Ramachandran Region Analysis:")
                        st.write(ramachandran_analysis)

        with tab2:
            st.subheader("Structural Comparison")
            uploaded_pdb1 = st.file_uploader("Upload the first PDB file", type=["pdb"], key="pdb1_compare")
            uploaded_pdb2 = st.file_uploader("Upload the second PDB file", type=["pdb"], key="pdb2_compare")
            pdb_data1_compare = None
            pdb_data2_compare = None
            if uploaded_pdb1 is not None and uploaded_pdb2 is not None:
                pdb_data1_compare = uploaded_pdb1.read().decode("utf-8")
                pdb_data2_compare = uploaded_pdb2.read().decode("utf-8")
                aligned_structure, rmsd = superimpose_structures(pdb_data1_compare, pdb_data2_compare)
                if aligned_structure:
                    st.success(f"Structures superimposed! RMSD: {rmsd:.2f} Å")
                    colA, colB = st.columns(2)
                    with colA:
                        st.markdown("**Structure 1**")
                        show_3d_structure(
                            pdb_data1_compare,
                            style=controls['render_style'],
                            highlight_ligands=controls['show_ligands'],
                            width=400,
                            height=400
                        )
                    with colB:
                        st.markdown("**Aligned Structure 2**")
                        show_3d_structure(
                            aligned_structure,
                            style=controls['render_style'],
                            highlight_ligands=controls['show_ligands'],
                            width=400,
                            height=400
                        )
                else:
                    st.error("Superimposition failed.")

if __name__ == "__main__":
    main()
