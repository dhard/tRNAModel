import argparse
import numpy as np
from random import Random

from evolver import Evolver
from models import Bit_Model, Test_Bit_Model, Fifty_Bit_Model, Ring_Bit_Model

if __name__ == "__main__":
    model_names = ["bit", "test", "50bit", "Ring"]
    parser = argparse.ArgumentParser(description="Evolve tRNAs and AARSs in a cmc model.")

    parser.add_argument("-g", "--population-size", type=int, default=100,
                        help="The size of the population to simulate (default: 100).")

    parser.add_argument("-aa", "--amino-acids", type=int, default=8,
                        help="The number of amino acids to model (default: 8).")

    parser.add_argument("-s", "--site-types", type=int, default=10,
                        help="The number of site types to model (default: 10).")

    parser.add_argument("-p", "--phi", type=float, default=0.5,
                        help="The intensity of selection on messages (default: 0.5).")

    parser.add_argument("-mm", "--message-mu", type=float, default=.01,
                        help="The rate of mutation for the messages (default: .01).")

    parser.add_argument("-m", "--mu", type=float, default=1e-4,
                        help="The rate of mutation for the SSWM model (default: 1e-4).")

    parser.add_argument("-w", "--site_weights", type=str, default=None,
                        help="The weight of each site represented as a string (e.g. 11111 for "
                        "a model with 5 sites) (defaults to all ones)")

    parser.add_argument("--seed", type=int, default=None,
                        help="The seed to use in the random number generator, defaults to a random value.")

    parser.add_argument("-n", "--model-name", choices=model_names, default=model_names[0],
                        help="The model type to be ran.")

    # For bit model
    parser.add_argument("-l", "--trna-length", type=int, default=3, 
                        help="The length of tRNAs in the bit model (default: 3).")
    parser.add_argument("-r", "--aars-length", type=int, default=3,
                        help="The length of AARSs in the bit model (default: 3).")
    parser.add_argument("-t", "--trnas", type=int, default=5, 
                        help="The number of tRNAs to include in the gnetic code (default: 5).")
    parser.add_argument("-a", "--aarss", type=int, default=5,
                        help="The number of AARSs to include in the gnetic code (default: 5).")


    args = parser.parse_args()

    rng = Random()
    rng.seed(args.seed)

    model = None

    if args.model_name == "bit":
        model = Bit_Model(args, rng)
    elif args.model_name == "test":
        model = Test_Bit_Model(args, rng)
    elif args.model_name == "50bit":
        model = Fifty_Bit_Model(args, rng)
    elif args.model_name == "Ring":
        model = Ring_Bit_Model(args, rng)

    evolver = Evolver(model.get_initial_code(), model.get_site_types(),
                      model.get_message_mutation_matrix(), args.population_size, rng)

    steps = 0
    mutations = 0
    last_code = model.get_initial_code()
    while mutations < 10:
        if evolver.frozen:
            print "Code is frozen, stopping..."
            break

        if last_code != evolver.current_code:
            last_code = evolver.current_code
            mutations += 1

            print steps, last_code

            #print last_code.codon_trna_mapping
            #print last_code.trna_aars_mapping
            #print last_code.aars_aa_mapping

            if 0 not in last_code.trnas and len(set(last_code.trnas)) == len(last_code.trnas):
                quit()

        evolver.step_time()
        steps += 1
