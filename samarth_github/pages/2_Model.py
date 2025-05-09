import streamlit as st
import requests
from Bio.PDB import PDBParser, PPBuilder
from io import StringIO
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import py3Dmol

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

def predict_structure_with_esmfold(sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post(url, headers=headers, data=sequence)
    if response.status_code == 200:
        return response.text
    else:
        st.error(f"ESMFold prediction failed: {response.text}")
        return None

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

def show_3d_structure(pdb_data, style='cartoon', highlight_ligands=True):
    view = py3Dmol.view(width=800, height=500)
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
    st.components.v1.html(view._make_html(), height=500, width=800)

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

def get_plddt_from_pdb(pdb_string):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pred", StringIO(pdb_string))
    b_factors = [atom.get_bfactor() for atom in structure.get_atoms()]
    if b_factors:
        return sum(b_factors)/len(b_factors)
    else:
        return None

# ----------------------
# UI Components
# ----------------------
def sidebar_controls():
    with st.sidebar:
        st.image("https://media.istockphoto.com/id/1390037416/photo/chain-of-amino-acid-or-bio-molecules-called-protein-3d-illustration.jpg?s=612x612&w=0&k=20&c=xSkGolb7TDjqibvINrQYJ_rqrh4RIIzKIj3iMj4bZqI=", width=400)
        st.title("Protein Molecule Mosaic")
        input_mode = st.radio(
            "Input Mode:",
            ["PDB ID", "Upload PDB", "Predict from Sequence (ESMFold)"],
            help="Choose how to provide a protein structure"
        )
        render_style = st.selectbox(
            "Rendering Style:",
            ["cartoon", "surface", "sphere"],
            index=0,
            help="Choose molecular representation style"
        )
        st.markdown("---")
        st.markdown("**Ligand Display Options**")
        show_ligands = st.checkbox("Highlight Ligands", True)
        return {
            'input_mode': input_mode,
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
    controls = sidebar_controls()
    col1, col2 = st.columns([3, 1])

    with col1:
        st.header("Protein Palette")
        st.markdown("**Load a protein structure:**")
        pdb_data = None
        source = None

        # --- NEW: Handle input modes ---
        if controls['input_mode'] == "PDB ID":
            pdb_id = st.text_input("Enter PDB ID (optional):").upper()
            if pdb_id:
                pdb_data = fetch_pdb_data(pdb_id)
                source = "pdbid"
                if pdb_data:
                    st.success(f"PDB ID {pdb_id} loaded from RCSB.")
        elif controls['input_mode'] == "Upload PDB":
            uploaded_pdb = st.file_uploader("Upload a PDB file", type=["pdb"])
            if uploaded_pdb is not None:
                pdb_data = uploaded_pdb.read().decode("utf-8")
                source = "upload"
                st.success("PDB file uploaded and loaded.")
        elif controls['input_mode'] == "Predict from Sequence (ESMFold)":
            example_seq = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
            seq = st.text_area("Input protein sequence (single-letter, no spaces):", value=example_seq, height=150)
            if st.button("Predict Structure with ESMFold"):
                with st.spinner("Predicting structure..."):
                    pdb_data = predict_structure_with_esmfold(seq)
                    source = "esmfold"
                    if pdb_data:
                        st.success("Structure predicted with ESMFold.")

        # --- Continue as before ---
        if pdb_data:
            st.subheader("3D Structure Viewer")
            show_3d_structure(
                pdb_data,
                style=controls['render_style'],
                highlight_ligands=controls['show_ligands']
            )
            if source == "esmfold":
                plddt = get_plddt_from_pdb(pdb_data)
                if plddt:
                    st.info(f"ESMFold average pLDDT confidence: {plddt:.2f}")

            with st.expander("Ramachandran Plot"):
                phi_psi = get_phi_psi_angles(pdb_data)
                if phi_psi:
                    fig = plot_ramachandran(phi_psi)
                    st.pyplot(fig)
                    stats = ramachandran_region_analysis(phi_psi)
                    st.markdown(f"""
                    **Ramachandran Plot Analysis**
                    - **Total residues:** {stats['total']}
                    - **Favored region:** {stats['favored']:.1f}%
                    - **Allowed region:** {stats['allowed']:.1f}%
                    - **Outlier region:** {stats['outlier']:.1f}%
                    """)
                else:
                    st.warning("Unable to generate Ramachandran plot. Please check the PDB input.")

            # --- Mutation Simulator ---
            with st.expander("🧬 Mutation Simulator (In Silico)"):
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("temp", StringIO(pdb_data))
                residues = [
                    (res.parent.id, res.id[1], res.get_resname())
                    for res in structure.get_residues()
                    if res.id[0] == ' '
                ]
                if residues:
                    residue_options = [
                        f"Chain {chain} Residue {resnum} ({resname})"
                        for chain, resnum, resname in residues
                    ]
                    selected = st.selectbox("Select residue to mutate:", residue_options)
                    selected_idx = residue_options.index(selected)
                    sel_chain, sel_resnum, sel_resname = residues[selected_idx]

                    aa_list = [
                        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                        "THR", "TRP", "TYR", "VAL"
                    ]
                    new_resname = st.selectbox("Mutate to:", aa_list, index=aa_list.index(sel_resname) if sel_resname in aa_list else 0)
                    if st.button("Mutate and Analyze"):
                        mutated_pdb = mutate_residue_in_pdb(pdb_data, sel_chain, sel_resnum, new_resname)
                        st.success(f"Residue mutated: Chain {sel_chain} {sel_resname}{sel_resnum} → {new_resname}{sel_resnum}")

                        st.subheader("Mutated Structure")
                        show_3d_structure(mutated_pdb, style=controls['render_style'], highlight_ligands=controls['show_ligands'])

                        st.subheader("Ramachandran Plot (Mutated)")
                        phi_psi_mut = get_phi_psi_angles(mutated_pdb)
                        if phi_psi_mut:
                            fig_mut = plot_ramachandran(phi_psi_mut)
                            st.pyplot(fig_mut)
                            stats_mut = ramachandran_region_analysis(phi_psi_mut)
                            st.markdown(f"""
                            **Ramachandran Plot Analysis (Mutated)**
                            - **Total residues:** {stats_mut['total']}
                            - **Favored region:** {stats_mut['favored']:.1f}%
                            - **Allowed region:** {stats_mut['allowed']:.1f}%
                            - **Outlier region:** {stats_mut['outlier']:.1f}%
                            """)
                        else:
                            st.warning("Unable to generate Ramachandran plot for mutated structure.")
                else:
                    st.info("No residues found for mutation.")

    with col2:
        st.header("Protein Dynamics")
        if pdb_data:
            with st.expander("Ligand Information"):
                ligands = extract_ligands(pdb_data)
                st.write(f"**Ions:** {len(ligands['ion'])}")
                st.write(f"**Ion Names:** {', '.join(ligands['ion'])}")
                st.write(f"**Monodentate Ligands:** {len(ligands['monodentate'])}")
                st.write(f"**Polydentate Ligands:** {len(ligands['polydentate'])}")
            with st.expander("Active Sites"):
                active_sites = predict_active_sites(pdb_data)
                st.write(f"**Predicted Active Sites ({len(active_sites)} residues):**")
                for site in active_sites:
                    st.write(f"{site['resname']} Chain {site['chain']} Residue {site['resnum']}")
                st.info("Active sites are predicted based on common catalytic residues (HIS, ASP, GLU, SER, CYS, LYS, TYR, ARG).")
            with st.expander("Ligand Type Visualization"):
                fig = visualize_ligand_counts(ligands)
                st.plotly_chart(fig)

if __name__ == "__main__":
    main()
