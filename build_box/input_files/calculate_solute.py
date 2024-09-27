import sys

if __name__ == "__main__":
    MOLECULE_NAME = 114.23 # MW MOLECULE_NAME
    water = 18.02 # MW water

    if len(sys.argv) != 3:
        print("Usage: python script.py <target_weight_percent> <n_water>")
        sys.exit(1)

    target = float(sys.argv[1])
    n_water = int(sys.argv[2])
    
    # calculate number of MOLECULE_NAME molecules
    # target = ( n_MOLECULE_NAME * MOLECULE_NAME ) / (( n_MOLECULE_NAME * MOLECULE_NAME ) + ( n_water * water))
    # Solve for n_MOLECULE_NAME:
    n_MOLECULE_NAME = int((target * n_water * water) / (MOLECULE_NAME * (1 - target)))
    
    print(f"The number of MOLECULE_NAME molecules is {n_MOLECULE_NAME}")