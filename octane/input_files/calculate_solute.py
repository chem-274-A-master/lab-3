import sys

if __name__ == "__main__":
    octane = 114.23 # MW octane
    water = 18.02 # MW water

    if len(sys.argv) != 3:
        print("Usage: python script.py <target_weight_percent> <n_water>")
        sys.exit(1)

    target = float(sys.argv[1])
    n_water = int(sys.argv[2])
    
    # calculate number of octane molecules
    # target = ( n_octane * octane ) / (( n_octane * octane ) + ( n_water * water))
    # Solve for n_octane:
    n_octane = int((target * n_water * water) / (octane * (1 - target)))
    
    print(f"The number of octane molecules is {n_octane}")