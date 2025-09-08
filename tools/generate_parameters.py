import csv
from argparse import ArgumentParser

def generate_parameters(filename):
    print(filename)

    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        data = [row for row in reader]

    name = data[0][1]
    singlet_energy_lowest = float(data[1][1])
    triplet_energy_lowest = float(data[2][1])
    kinetic_energy_factor = float(data[3][1])
    kinetic_energy_power = float(data[4][1])
    ionize_energy = []
    for row in data[5:]:
      if row[0] == '' and row[1] == '':
        break
      ionize_energy.append(float(row[1]))
    norbital = len(ionize_energy)

    diff_singlet_energy = ionize_energy[0]-singlet_energy_lowest
    diff_triplet_energy = ionize_energy[0]-triplet_energy_lowest
    singlet_energy = [e - diff_singlet_energy for e in ionize_energy[:norbital]]
    triplet_energy = [e - diff_triplet_energy for e in ionize_energy[:norbital]]
    kinetic_energy = [kinetic_energy_factor * (e ** kinetic_energy_power) for e in ionize_energy[:norbital]]

    outfile = name.lower() + "_nml.in"
    with open(outfile,"w") as fp:
        fp.write(f"!\n")
        fp.write(f"! Auto-generated file: Do not edit!\n")
        fp.write(f"!\n")
        fp.write(f"\n")
        fp.write(f"&param_orbital\n")
        fp.write(f" norbital = {norbital:0d}\n")
        fp.write(f" name = '{name}'\n")
        fp.write(f"/\n")
        fp.write(f"\n")
        fp.write(f"&parameters_per_orbitals\n")
        fp.write(f" {'energy_ionize':16s} = ")
        for io in range(norbital):
            fp.write(f" {ionize_energy[io]:8.2f}")
        fp.write(f"\n")
        fp.write(f" {'energy_singlet':16s} = ")
        for io in range(norbital):
            fp.write(f" {singlet_energy[io]:8.2f}")
        fp.write(f"\n")
        fp.write(f" {'energy_triplet':16s} = ")
        for io in range(norbital):
            fp.write(f" {triplet_energy[io]:8.2f}")
        fp.write(f"\n")
        fp.write(f" {'energy_kinetic':16s} = ")
        for io in range(norbital):
            fp.write(f" {kinetic_energy[io]:8.2f}")
        fp.write(f"\n")
        fp.write(f" {'number_electrons':16s} = ")
        for io in range(norbital):
            fp.write(f" {2:8d}")
        fp.write(f"\n")
        fp.write(f"/\n")

if __name__ == '__main__':
    # parse arguments
    parser = ArgumentParser()
    parser.add_argument('filename')

    args = parser.parse_args()

    generate_parameters(args.filename)
