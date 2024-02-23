import numpy as np


def viterbi(observations, states, start_prob, trans_prob, emm_prob):
    T = len(observations)
    N = len(states)

    viterbi_matrix = np.zeros((N, T))
    backpointer = np.zeros((N, T), dtype=int)

    # Initialization step
    for s in range(N):
        viterbi_matrix[s, 0] = start_prob[s] * emm_prob[s][observations[0]]
        backpointer[s, 0] = 0

    print("after Initialization step")
    print(viterbi_matrix)
    print(backpointer)

    # Recursion step
    for t in range(1, T):
        for s in range(N):
            viterbi_matrix[s, t] = max(
                viterbi_matrix[s_prime, t - 1] * trans_prob[s_prime][s] * emm_prob[s][observations[t]] for s_prime in
                range(N))
            backpointer[s, t] = np.argmax(
                [viterbi_matrix[s_prime, t - 1] * trans_prob[s_prime][s] for s_prime in range(N)])

    print("after Recursion step")
    print(viterbi_matrix)
    print("backpointer\n",backpointer)

    # Termination step
    best_path_prob = max(viterbi_matrix[s, T - 1] for s in range(N))
    best_path_pointer = np.argmax([viterbi_matrix[s, T - 1] for s in range(N)])

    print("after Termination step")
    print("best_path_prob\n",best_path_prob)
    print("best_path_pointer\n",best_path_pointer)

    # Path backtracking
    best_path = [best_path_pointer]
    for t in range(T - 1, 0, -1):
        best_path_pointer = backpointer[best_path_pointer, t]
        best_path.insert(0, best_path_pointer)

    print("after Path backtracking")
    print("best_path\n", best_path)

    # Convert state indices to state names
    best_path_states = [states[index] for index in best_path]

    return best_path_states, best_path_prob

# Example usage:
# Define your observations, states, start probabilities, transition probabilities, and emission probabilities here.

# Mapping the nucleotide sequence to indices
nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
sequence = 'CGAAAAAATCG'
observations = [nucleotide_to_index[nuc] for nuc in sequence]

# States
states = ['non_coding', 'coding']
state_to_index = {'coding': 0, 'non_coding': 1}


# Start probabilities (example values)
start_prob = [1, 0]  # Assuming more probability to start in non coding state


# Transition probabilities (example values)
trans_prob = [[0.8, 0.2],  # Stay in 'non_coding' or switch to 'coding'
              [0.4, 0.6]]  # Switch to 'non_coding' or stay in 'coding'

# Emission probabilities (example values)
emm_prob = [[0.2, 0.3, 0.3, 0.2], # Probabilities for 'A', 'C', 'G', 'T' in 'non_coding'
            [0.4, 0.2, 0.2, 0.2]]  # Probabilities for 'A', 'C', 'G', 'T' in 'coding'

best_path, best_path_prob = viterbi(observations, states, start_prob, trans_prob, emm_prob)
print("Best Path: ", best_path)
print("Best Path Probability: ", best_path_prob)
