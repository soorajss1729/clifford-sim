import streamlit as st
import numpy as np
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt

#tableau = np.full((2 * k, 2 * k + 1), np.nan, dtype=object)
#gates=[]
#-----------------------------
def anticom_rows_gen(k):
	found=False
	while not found:
		r1 = np.random.randint(0, 2, size=(2 * k + 1))
		r2= np.random.randint(0, 2, size=(2 * k + 1))
		symplectic_ip=np.sum(r1[:2*k:2] * r2[1:2*k:2] + r1[1:2*k:2] * r2[:2*k:2]) % 2 # [start : stop : step]
		if symplectic_ip==1:
			found=True
	return r1,r2
#-----------------------------	 
#tableau[0], tableau[1] = anticom_rows_gen(k)
#tableau[0]=[1,0,1,1,1,1,1,0,0]
#tableau[1]=[1,1,1,1,1,1,1,0,0]

#tableau[0]=[0,0,0,1,0,0,0]
#tableau[1]=[1,1,1,1,0,0,0]


#print('stp0:',tableau)
#-----------------------------
def extract_sweep_tabs(tableau, k,pair_index):
	sweep_tab_x = np.zeros((2, k), dtype=int)
	sweep_tab_z = np.zeros((2, k), dtype=int)
	sweep_tab_x[0, :k] = tableau[2*pair_index, 2*pair_index:2*pair_index+2*k:2]  # First k elements from even columns
	sweep_tab_z[0, :k] = tableau[2*pair_index, 2*pair_index+1:2*pair_index+2*k:2]  # Next k elements from odd columns
	sweep_tab_x[1, :k] = tableau[2*pair_index+1, 2*pair_index:2*pair_index+2*k:2]  # First k elements from even columns
	sweep_tab_z[1, :k] = tableau[2*pair_index+1, 2*pair_index+1:2*pair_index+2*k:2]  # Next k elements from odd columns
	print('sweep_x:',sweep_tab_x)
	print('sweep_z:',sweep_tab_z)
	print('extract_sweep_tabs done')
	return sweep_tab_x,sweep_tab_z
#-----------------------------

#-----------------------------
def clear_z_block(row, k, sweep_tab_x, sweep_tab_z, gates):
    for j in range(k):
        x_a = sweep_tab_x[row, j]  # Use from sweep_tab_x
        z_a = sweep_tab_z[row, j]  # Use from sweep_tab_z
        if z_a == 1:
            if x_a == 0:
                gates.append(("h", j))  # x_a <--> z_a
                sweep_tab_x[0, j], sweep_tab_z[0, j] = sweep_tab_z[0, j], sweep_tab_x[0, j]
                sweep_tab_x[1, j], sweep_tab_z[1, j] = sweep_tab_z[1, j], sweep_tab_x[1, j]
            else:
                gates.append(("s", j))  # z_a = x_a + z_a
                sweep_tab_z[0, j] ^= sweep_tab_x[0, j]
                sweep_tab_z[1, j] ^= sweep_tab_x[1, j]                
    print('sweep_x:',sweep_tab_x)
    print('sweep_z:',sweep_tab_z)
    print('gates:',gates)
#    return sweep_tab_x,sweep_tab_z
#-----------------------------
def sweep_x_to_pivot(row, k, sweep_tab_x, sweep_tab_z, gates):
     J = [j for j in range(k) if sweep_tab_x[row, j] == 1]
     print('J:',J)
     while len(J)>1:
          J = sorted(J)
          J_1=[]
          for i in range(0,len(J)-1,2):
               print(i)
               c, t = J[i], J[i+1]
               sweep_tab_x[0, t] ^= sweep_tab_x[0, c]
               sweep_tab_z[0, c] ^= sweep_tab_z[0, t]
               sweep_tab_x[1, t] ^= sweep_tab_x[1, c]
               sweep_tab_z[1, c] ^= sweep_tab_z[1, t]
               gates.append(("cx",c,t))
               J_1.append(c)
          if len(J) % 2 == 1:
               J_1.append(J[-1])
          J=J_1
     print('JJ:',J)
     print('sweep_x:',sweep_tab_x)
     print('sweep_z:',sweep_tab_z)
     print('gates:',gates)
     	
     if J[0]!=0:
          for (c,t) in [(0,J[0]),(J[0],0),(0,J[0])]:
               sweep_tab_x[0, t] ^= sweep_tab_x[0, c]
               sweep_tab_z[0, c] ^= sweep_tab_z[0, t]
               sweep_tab_x[1, t] ^= sweep_tab_x[1, c]
               sweep_tab_z[1, c] ^= sweep_tab_z[1, t]
               gates.append(("cx",c,t))
     print('sweep_x:',sweep_tab_x)
     print('sweep_z:',sweep_tab_z)
     print('gates:',gates)
#     return sweep_tab_x,sweep_tab_z
#-----------------------------
def sweep_second_row(sweep_tab_x, sweep_tab_z, k, gates):
     gates.append(("h",0))
     sweep_tab_x[0,0], sweep_tab_z[0,0] = sweep_tab_z[0,0], sweep_tab_x[0,0]
     sweep_tab_x[1,0], sweep_tab_z[1,0] = sweep_tab_z[1,0], sweep_tab_x[1,0]
     print('sweep_x:',sweep_tab_x)
     print('sweep_z:',sweep_tab_z)
     clear_z_block(1, k, sweep_tab_x, sweep_tab_z, gates)
     sweep_x_to_pivot(1, k, sweep_tab_x, sweep_tab_z, gates)
     gates.append(("h", 0))
     sweep_tab_x[0, 0], sweep_tab_z[0, 0] = sweep_tab_z[0, 0], sweep_tab_x[0, 0]
     sweep_tab_x[1, 0], sweep_tab_z[1, 0] = sweep_tab_z[1, 0], sweep_tab_x[1, 0]
#-----------------------------
def sweep_pair(tableau, k, step, gates):
	sweep_tab_x,sweep_tab_z=extract_sweep_tabs(tableau, k, step)
	clear_z_block(0, k, sweep_tab_x, sweep_tab_z, gates)
	sweep_x_to_pivot(0, k, sweep_tab_x, sweep_tab_z, gates)
	if not ((sweep_tab_x[1].sum()==0) and (sweep_tab_z[1].sum()==1) and (sweep_tab_z[1,0] == 1)):
		print('Step 4 Executing!!!')
		sweep_second_row(sweep_tab_x, sweep_tab_z, k, gates)
	print('!! Iteration End !!')
	return sweep_tab_x,sweep_tab_z
#sweep_tab_x,sweep_tab_z=sweep_pair(tableau, k, 0, gates)
#print('sweep_x:',sweep_tab_x)
#print('sweep_z:',sweep_tab_z)

def random_clifford(k):
	tableau = np.full((2 * k, 2 * k + 1), np.nan, dtype=object)
	gates = []
	remaining = k
	for i in range(k):
		print('### Iteration ### :',i)
		r1, r2 = anticom_rows_gen(remaining)
		start_row = 2*i
		start_col = 2*i
		tableau[start_row,     start_col : start_col + 2*remaining+1] = r1 # start:stop:step=1
		tableau[start_row + 1, start_col : start_col + 2*remaining+1] = r2		
		print('hooooo',tableau)
		sweep_pair(tableau, remaining, i, gates)
		remaining -= 1
	   # Create the reversed and daggered gate list
	daggered_gates = []
	for gate in reversed(gates):
		gate_type, *params = gate
		if gate_type == "h":
			daggered_gates.append(("h", *params))  # H is self-inverse
		elif gate_type == "s":
			daggered_gates.append(("sd", *params))  # Daggered S (sd)
		elif gate_type == "cx":
			daggered_gates.append(("cx", *params))  # CX is self-inverse		
	return gates, daggered_gates, tableau


def create_clifford_circuit(k, daggered_gates):
    qc = QuantumCircuit(k)
    
    for gate in daggered_gates:
        gate_type = gate[0]
        
        if gate_type == "h":
            qc.h(gate[1])
        elif gate_type == "sd":
            qc.sdg(gate[1])
        elif gate_type == "cx":
            qc.cx(gate[1], gate[2])
        else:
            raise ValueError(f"Unknown gate type: {gate_type}")
    return qc

# Streamlit UI
st.title("Interactive Random Clifford Operator Simulator")
st.write("Generate and visualize a random Clifford operator.")

# User input: Number of qubits (k)
k = st.slider("Select the number of qubits (k):", min_value=1, max_value=10, value=5)

if st.button("Generate Random Clifford Operator"):
    # Generate Clifford Gates and Tableau
    gates_sequence, daggered_gates_sequence, final_tableau = random_clifford(k)
    
    st.subheader("Generated Clifford Operator Gates (Compact List)")
    formatted_gates = ", ".join([str(gate) for gate in daggered_gates_sequence])
    st.text(formatted_gates)
    
#    st.subheader("Stabilizer Tableau")
#    st.write(final_tableau)
    
    # Generate and display Clifford Circuit
    daggered_clifford_circuit = create_clifford_circuit(k, daggered_gates_sequence)
    st.subheader("Random Clifford Circuit")
    st.pyplot(daggered_clifford_circuit.draw(output='mpl'))




